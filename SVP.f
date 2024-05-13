!to add in the run_star_extras.f

   module run_star_extras      
    use star_lib
    use star_def
    use const_def
    use math_lib
    use chem_def
    use chem_lib
    use ionization_lib
    use MOD_SVP, only: initialize_SVP, g_rad_SVP
    implicit none
    
    
    
    

    subroutine SVP_routine( id, &
            nz, nzlo, nzhi, nc, m, kmax_rad_accel, X, A, &
            class_chem_id, net_iso, op_mono_factors, &
            L_face, rho_face, r_face, T_face, alfa_face, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening, log10_g_rad, g_rad, &
            rad_accel_face, ierr)
         
         
         ! inputes that MESA uses in its radiative accelerations, and are considered outside the SVP routine
         integer, intent(in) :: id
         integer, intent(in) :: nz, nzlo, nzhi, nc, m, kmax_rad_accel, class_chem_id(:), net_iso(:), &
                                min_Z_for_radaccel, max_Z_for_radaccel
         real(dp), intent(in) :: min_T_for_radaccel, max_T_for_radaccel
         real(dp), dimension(:), intent(in) :: A, L_face, rho_face, r_face, T_face, alfa_face, op_mono_factors
         ! A ->atomic mass number; L_face -> local luminosaty; rho_face -> local density  
         ! T_face -> local Temperature; alfa_face -> interpolation parameter; op_mono_factors -> OP mono opacyties factores
         ! e.g., T_face(k) = alfa_face(k)*T(k) + (1d0-alfa_face(k))*T(k-1)
         real(dp), dimension(:,:), intent(in) :: X
         logical, intent(in) :: screening
         ! radiative accelerations
         real(dp), dimension(:,:), intent(out) :: log10_g_rad, g_rad, rad_accel_face
         integer, intent(out) :: ierr
         
         
         ! total mass and radius, effective Temperature, and initial hidrogen abundance
         real(dp) :: mtot, rtot, teff, Xi
         ! vg_rad_svp and dvg_rad_svp parameters used to compute g_rad and its derevative
         ! ychim is the chemical composition per mole
         real(dp), dimension(100) :: vg_rad_svp, ychim
         real(dp), dimension(100,200) :: dvg_rad_svp
         ! total number of elements followed by atomic diffusion.
         integer :: nchim
         integer  :: i,j,k
         !Number of protons (atomic number), number of neutrons, and number of particles in the nucleo
         integer, dimension(100) :: Zi,Ni,Z_Ni
         !Atomic weight
         real(dp), dimension(100) :: Aj
         !Mean charge of each particle at each zone
         real(dp), dimension(100,4000) :: typical_charge
         !id number of the element in MESA
         integer, dimension(100)  :: elem_id
         !element name (e.g. C12 is carbon 12)
         character(8), dimension(100):: elem_name
         !electron density 
         real(dp) :: Mnel, in_mu_e
         !mass, radius, temperature and local gravity at that region 
         real(dp) :: mass,ray, temp, grav
         !corrent version of svp
         character(64) :: version_svp
         
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
               
         
         
         !Set the outpots deffault as zero 
         vg_rad_svp=0.0
         dvg_rad_svp=0.0
         g_rad(:,:) = 0.0
         log10_g_rad(:,:) = -99
         rad_accel_face(:,:) = 0.0
         !SVP only valid in the MS (between ZAMS and TAMS)
         Xi=1-s% initial_z-s% initial_y
         if (s%X(nz)<=1d-5)  then
            return
         else if (abs((Xi-s%center_h1)/Xi)<=0.01) then
            return
         end if
         !defining internal variavais (to avoid NAN or random values)
		     elem_id(:)=-1  
         elem_name(:)= ''
         typical_charge(:,:)=0
         Aj(:)=0;
         Zi(:)=0;Ni(:)=0;Z_Ni(:)=0
         ychim(:)=0.0
         Mnel=0
         !defining number of elements follows in atomic diffusion, total stellar mass, radius, luminosaty and effective temperature
		     nchim = size(class_chem_id)
         mtot = s% mstar	
         rtot = s% photosphere_r*Rsun
         teff  = s% Teff	 
         !SVP only works in stellar masses betwen 1.0 Msun and 10.0 Msun
         if (mtot<1.0 .and. mtot>10.0) then
            return
         end if 
         !-----------------------------------------------------------------------------------------------------------------------------
         !to obtain the id number of the element in MESA
         i=1
         j=1
         do k=1, size(s% net_iso)
           if  (s% net_iso(k)>0 .and. k==class_chem_id(i)) then
               elem_id(i)=j
               i=i+1
               j=j+1
           else if (s% net_iso(k)>0) then
               j=j+1
           end if
           if (i>nchim) then
              exit
           end if
         end do
         !to obtain the Atomic weight, the Atomic number, number of protons and Atomic mass number
         do k=1, nchim
            Aj(k)=chem_isos% W(class_chem_id(k))
            Zi(k)=chem_isos% Z(class_chem_id(k))
            Ni(k)=chem_isos% N(class_chem_id(k))
            Z_Ni(k)=Zi(k)+Ni(k)
            elem_name(k)= chem_isos% name(class_chem_id(k))          
         end do
         !To obtain the mean charge of each element on all the stellar profile
         do i=1, nz
            do j=1, nchim
               typical_charge(j,i) = eval_typical_charge( class_chem_id(j), s% abar(i), exp(s% lnfree_e(i)), &
                                                          s% T(i), s% lnT(i)/log(10.), s% rho(i), s% lnd(i)/log(10.))
            end do
         end do
         !-----------------------------------------------------------------------------------------------------------------------------
         !initialize the SVP needs to be called at the ZAMS and only one time per model
         !any of the s% x_integer_ctrl(k) can be use by deffault they are set to zero.
         if (s% x_integer_ctrl(4)==0) then
            version_svp='datai_SVP_v1'
            call initialize_SVP(mtot,nchim,elem_name(1:nchim),Aj(1:nchim),Zi(1:nchim),version_svp)
            s% x_integer_ctrl(4)=-1
         end if
         !-----------------------------------------------------------------------------------------------------------------------------
         do i=1, nz 
         		!Calculation of electron density
            if (s% x_logical_ctrl(7)) then   
               !full ionisation aproximation 
               Mnel= s%rho(i) * avo*(1+s% X(i))/2
            else
               !Partial ionisation consideration 
               Mnel= s%rho(i) * avo*exp(s%lnfree_e(i))
            end if
            !profile parameters: mass, radius, temperature, and local gravity
            mass = 1.0 - s% q(i)   !mass above the layer (0 at the surface, 1 at the center)
            ray = s% r(i)
            temp = s% T(i)
            grav = s% grav(i)
            !chemical composition per mole
            do j=1, nchim
               ychim(j) =s% xa(elem_id(j),i)/Aj(j)
            end do
         		!SVP computation of grad
            vg_rad_svp(1:nchim)=0.
            dvg_rad_svp(1:nchim,:)=0.
            !call for SVP routine provided by Alecian & LeBlanc (2020)
            call g_rad_SVP(nchim,elem_name(1:nchim),Aj(1:nchim),Zi(1:nchim),mtot,rtot,teff, &
                            mass,ray,temp,Mnel,ychim(1:nchim), &
                            vg_rad_svp(1:nchim),dvg_rad_svp(1:nchim,:))
                            
            g_rad(1:nchim,i)=vg_rad_svp(1:nchim)
            rad_accel_face(1:nchim,i)=vg_rad_svp(1:nchim)
            !to avoid NAN or inf in the log10(g_rad)
            do j=1, nchim
               if (vg_rad_svp(j)==0.0) then
                  log10_g_rad(j,i)=-99
               else
                  log10_g_rad(j,i)=log10(vg_rad_svp(j))
               end if
            end do
         end do
         
      end subroutine SVP_routine



  end module run_star_extras
