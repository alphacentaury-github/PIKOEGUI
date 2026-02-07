* In GUI, file unit numbers are fixed by default values.
* In GUI, shell commands are used to run external programs 
  at the input file locations. 
  In other words, "current working directory" are 
  (1) location of "cnt" input files when calling pikoe. 
      Thus, all paths are written w.r.t. the cnt file location. 
  (2) location of "omget.inp" file when calling omget. 
  Thus, GUI generated input "cnt" files and original "cnt" files
  can have different unit numbers and file paths.  

*default input/output file units in GUI 

input:
  cnt files : "./cnt/_test.cnt" 
  kibelm=11 : elemtary interaction file (for ielm=3,4) , '../elem/nnampFL.dat'
  ish = 12  : external s.p. wave file , './in_ish.dat'
  ipot0 = 20 : optical potetential file for particle 0  '../pot/EDAD1p12C_e.dat'
  ipot1 = 21 : optical potetential file for particle 1  '../pot/EDAD1p11B_e.dat'
  ipot2 = 22 : optical potetential file for particle 2  '../pot/EDAD1p11B_e.dat'

output: 
  kibbs = 13 (for kibbs=1) './out_kibbs.dat' 
  kibout = 6  : basic terminal output file '../outlist/12Cp2pTDXnorm.outlist' 
  kibtbl = 10 : tbl output file '../tbl/tbl_12Cp2pTDXnorm.dat' 
  kibtmd = 15 : './out_kibtmd.dat'
  kiblg = 16 :  './out_kiblg.dat'
  kibpx = 17 :  './out_kibpx.dat'
  kibtr = 18 :  './out_kibtr.dat'
  kibtl = 19 :  './out_kibtl.dat'
  ivar  = 14 :input?output? for kinematical output(TDX,QDX,...) files, './in_ivar.dat' 
        = 24
        = 34
        = 44            
