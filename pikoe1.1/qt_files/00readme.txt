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

* TO DO LIST
(1) make GUI to read "cnt" file directly.
(2) make omget to prepare direct input potential file for pikoe.  

#---------------------
When file is required 
#---------------------
(1) ish > 9   <--> ish = 12 in gui and filename 
(2) ielm=3 or 4 <---> kibelm = 11 in gui 
(3) ipot0 > 9 <--> ipot0 = 20 in gui 
(4) ipot1 > 9 <--> ipot1 = 21 in gui 
(5) ipot2 > 9 <--> ipot2 = 22 in gui 
(6) kibbs > 0 <--> kibbs = 13 in gui 
(7) kibout > 0 <--> kibout = 6 in gui 
(8) kibtbl >0 <--> kibtbl = 10 in gui   
(9) kibtmd >0 <--> kibtmd = 15 in gui 
(10) ivar=9 and kiblg > 0 <--> kiblg = 16 in gui 
(11) ivar=9 and kibpx > 0 <--> kibpx = 17 
(12) ivar=9 and kibtr > 0 <--> kibtr = 18  
(13) ivar=9 and kibrl > 0 <--> kibtl =19 
(14) 9 < ivar < 20  <--> ivar = 14 in gui 
(15) 19 < ivar < 30 <--> ivar = 24 in gui 
(16) 29 < ivar < 40 <--> ivar = 34 in gui 
(17) 39 < ivar      <--> ivar = 44 in gui 

