# SVP_routine-MESA_implementation
SVP routine for implementation in mesa-r12778 and mesa-24.08.1.

This is the SVP routine provided by Alecian & LeBlanc (2020 and 2024) in (https://gradsvp.obspm.fr/index.html) that I implemented in MESA (Modules for Experiments in Stellar Astrophysics, mesa-r12778). 

Instructions.txt - This is the instructions on how to implement SVP in the MESA code.

SVP.f - The routine that needs to be implemented in run_star_extra.f of MESA, after the modifications in the source code

SVP_2020.zip - it contains the SVP tables (version released in 2020), and the main routines of SVP that are called in SVP.f, the ''mod_svp.f'' and ''svp_codes.f''. The mod_svp.f should be called in the makefile that compiles before running the model. 

SVP_2024.zip - it contains the SVP tables (version released in 2024), and the main routines of SVP that are called in SVP.f, the ''mod_svp.f'' and ''svp_codes.f''. The mod_svp.f should be called in the makefile that compiles before running the model. (It updates the interpolation routine used)
               
Any questions you can send me an email to nammo@dtu.dk

