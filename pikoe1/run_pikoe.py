# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 11:32:14 2025

@author: User
"""
import sys
import os
import shutil
from io import StringIO
from subprocess import (call,Popen)
import subprocess
import numpy as np
from numpy import loadtxt

import json
import pickle

class PIKOE_input:
    float_variables = ['ZP', 'AP','ZA','AA','ELAB',
                       'EBIND','ZSP','ASP','BETASP',
                       'FJ','FL','SFAC',
                       'V0LS','RS','AS','FKNCUT',
                       'RC','RCL','A0C',
                       'VARMIN','VARMAX','DVAR',
                       'THXMIN','THXMAX','DTHX',
                       'PHXMIN','PHXMAX','DPHX',
                       'ET2MIN','ET2MAX','DET2',
                       'PH2MIN','PH2MAX','DPH2',
                       'RMAX','DR',
                       'FACV0', 'FACW0', 'FACVS0', 'FACWS0','BETA0', 
                       'FACV1', 'FACW1', 'FACVS1', 'FACWS1','BETA1', 
                       'FACV2', 'FACW2', 'FACVS2', 'FACWS2','BETA2', 
                       ]
    int_variables = ['LIMFS','IONS','IFRM','IMIR','ICAL',
                     'IKIN','ICTREIN','ISH','ICTRM','NOD','KIBBS',
                     'ICTRC','ICTRCL',
                     'IBMC','IBMS','ICTRS',
                     'LMAX0','LMAX1','LMAX2',
                     'IVAR','IEX','IXUNT','KUNT','IVVAR','IVTHX','IVPHX',
                     'IVET2','IVPH2',
                     'KIBTBL','KIBOUT','KIBTMD','KIBLG','KIBPX','KIBTR','KIBTL',
                     'IELM','KIBELM','IONSH','KINELM','IELMEDG',
                     'NGR','NGTH','NGPH','NGK1','NGPH1Q',
                     'IPOT0','IMS0','IEDG0', 'IPOT1','IMS1','IEDG1',
                     'IPOT2','IMS2','IEDG2'
                     ]
    
    str_variables = ['COMMENT' ]
    def __init__(self):
        #---default values 
        self.data = { 
            'HEADERS': [ [10,'unknown','./tbl_12Cp2pTDXnorm.dat'],
                         [11,'old','../elem/nnampFL.dat'],
                         [12,'old','../pot/EDAD1p12C_e.dat'],
                         [13,'old','../pot/EDAD1p11B_e.dat'],
                         [6,'unknown','./12Cp2pTDXnorm.outlist']],
            'COMMENT': '12C(p,2p)11B_gs@392MeV DWIA TDX normal',
            'LIMFS': 1000,'IONS' : 0,'IFRM': 0,'IMIR':0,'ICAL':1,
            'ZP': 1.0, 'AP':1.007825,'ZA':6.0,'AA':12.0,
            'IKIN':0,'ELAB':392.0,'ICTREIN': 0,
            'ISH':1,'EBIND':15.96,'ZSP':1.0,'ASP':1.007825,
            'BETASP':0.85, 'ICTRM': 1,
            'FJ': 1.5,'FL':1.0,'SFAC':1.77,'NOD': 0,'KIBBS':0, 
            'IBMC':0,'RC':1.35,'ICTRC':1,'A0C':0.65,'RCL':1.35,'ICTRCL':1,
            'IBMS':0,'V0LS':8.2,'RS':1.35,'ICTRS':1,'AS':0.65,
            'LMAX0':60,'LMAX1':60,'LMAX2':60,
            'IVAR':1,'IEX':0,'FKNCUT':2.0,'IXUNT':1,'KUNT':1,
            'IVVAR':0,'VARMIN':251.0,'VARMAX':255.0,'DVAR':10.0,
            'IVTHX':0,'THXMIN':32.5,'THXMAX':180.0,'DTHX':10.0, 
            'IVPHX':0,'PHXMIN':0.0,'PHXMAX':40.0,'DPHX':10.0,
            'IVET2':1,'ET2MIN':0.0,'ET2MAX':180.0,'DET2':0.5,
            'IVPH2':0,'PH2MIN':180.0,'PH2MAX':360.0,'DPH2':10.0,
            'KIBTBL':10,'KIBOUT':6,'KIBTMD':0,'KIBLG':0,'KIBPX':0,'KIBTR':0,'KIBTL':0,
            'IELM':4,'KIBELM':11,'IONSH':1,'KINELM':0,'IELMEDG':1,
            'RMAX':15.0,'DR':0.1,'NGR':30,'NGTH':30,'NGPH':40,'NGK1':0,'NGPH1Q':0,
            'IPOT0' : 12 ,'FACV0': 1.0, 'FACW0': 1.0 , 'FACVS0': 1.0 , 'FACWS0':1.0,
            'BETA0' : -0.85, 'IMS0': 0, 'IEDG0' : 1, 
            'IPOT1' : 13 ,'FACV1': 1.0, 'FACW1': 1.0 , 'FACVS1': 1.0 , 'FACWS1':1.0,
            'BETA1' : -0.85, 'IMS1': 0, 'IEDG1' : 1, 
            'IPOT2' : 13 ,'FACV2': 1.0, 'FACW2': 1.0 , 'FACVS2': 1.0 , 'FACWS2':1.0,
            'BETA2' : -0.85, 'IMS2': 0, 'IEDG2' : 1
                } 
        self.file_txt ='' 
        
    def set_data(self,**attributes):
        """update dictionary item values """
        for item,value in attributes.items():
            if item in self.str_variables: 
                    self.data[item]=value
            elif item in self.int_variables:
                    self.data[item]=int(float(value)) 
            elif item in self.float_variables:
                    self.data[item]=float(value)
            elif item in ['HEADERS','POTS']:
                    self.data[item]=value
            else :
                    print('invalid input item name %s\n'%item)
                    
    def write_txt(self,):
        """convert dictionary into pikoe input text files"""
        out_txt =' **** ppN control data **** \n'
        for  ff in self.data['HEADERS']:
            # "A1,I3,1X,A8,2X,A50" format 
            out_txt += ' {:3d}:{:8s}::{:50s}\n'.format(ff[0],ff[1],ff[2])
        out_txt += ' 999: \n'
        out_txt += '---- input ----\n'  
        out_txt += self.data['COMMENT']+'\n'  #L1
        out_txt += '{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}\n'.format(self.data['LIMFS'],self.data['IONS'],
                        self.data['IFRM'],self.data['IMIR'],self.data['ICAL']) #L2
        out_txt += '{:<5.2f}{:<10g}{:<5.2f}{:<10g}\n'.format(
            self.data['ZP'],self.data['AP'],self.data['ZA'],self.data['AA']) #L3
        out_txt += '{:<5d}{:<10g}{:<5d}\n'.format(
            self.data['IKIN'],self.data['ELAB'],self.data['ICTREIN']) #L4
        out_txt += '{:<5d}{:<10g}{:<5g}{:<10g}{:<10g}{:<5d}\n'.format(
            self.data['ISH'],self.data['EBIND'],
            self.data['ZSP'],self.data['ASP'],self.data['BETASP'],self.data['ICTRM']) #L5
        out_txt += '{:<5g}{:<5g}{:<10g}{:<5d}{:<5d}\n'.format(
            self.data['FJ'],self.data['FL'],self.data['SFAC'],self.data['NOD'],self.data['KIBBS']) #L6 
        out_txt += '{:<5d}{:<10g}{:<5d}{:<10g}{:<10g}{:<5d}\n'.format(
            self.data['IBMC'],self.data['RC'],self.data['ICTRC'],
            self.data['A0C'],self.data['RCL'],self.data['ICTRCL'] ) #L7
        out_txt += '{:<5d}{:<10g}{:<10g}{:<5d}{:<10g}\n'.format(
            self.data['IBMS'],self.data['V0LS'],self.data['RS'],self.data['ICTRS'],self.data['AS'] ) #L8
        out_txt += '{:<5d}{:<5d}{:<5d}\n'.format(self.data['LMAX0'],self.data['LMAX1'],self.data['LMAX2']) #L9
        out_txt += '{:<5d}{:<5d}{:<5g}{:<5d}{:<5d}\n'.format( 
            self.data['IVAR'],self.data['IEX'],self.data['FKNCUT'],self.data['IXUNT'],self.data['KUNT']) #L10 
        out_txt += '{:<5d}{:<10g}{:<10g}{:<10g}\n'.format( 
            self.data['IVVAR'],self.data['VARMIN'],self.data['VARMAX'],self.data['DVAR']) #L11 
        out_txt += '{:<5d}{:<10g}{:<10g}{:<10g}\n'.format( 
            self.data['IVTHX'],self.data['THXMIN'],self.data['THXMAX'],self.data['DTHX']) #L12 
        out_txt += '{:<5d}{:<10g}{:<10g}{:<10g}\n'.format( 
            self.data['IVPHX'],self.data['PHXMIN'],self.data['PHXMAX'],self.data['DPHX']) #L13
        out_txt += '{:<5d}{:<10g}{:<10g}{:<10g}\n'.format( 
            self.data['IVET2'],self.data['ET2MIN'],self.data['ET2MAX'],self.data['DET2']) #L14 
        out_txt += '{:<5d}{:<10g}{:<10g}{:<10g}\n'.format( 
            self.data['IVPH2'],self.data['PH2MIN'],self.data['PH2MAX'],self.data['DPH2']) #L15 
        out_txt += '{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}\n'.format(
            self.data['KIBTBL'],self.data['KIBOUT'],self.data['KIBTMD'],
            self.data['KIBLG'],self.data['KIBPX'],self.data['KIBTR'],self.data['KIBTL']) #L16 
        out_txt += '{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}\n'.format(
            self.data['IELM'],self.data['KIBELM'],self.data['IONSH'],self.data['KINELM'],self.data['IELMEDG']) #L17
        out_txt += '{:<10g}{:<10g}{:<5d}{:<5d}{:<5d}{:<5d}{:<5d}\n'.format(
            self.data['RMAX'],self.data['DR'],self.data['NGR'],self.data['NGTH'],self.data['NGPH'],
            self.data['NGK1'],self.data['NGPH1Q']) #L18
        out_txt += '{:<5d}{:<5g}{:<5g}{:<5g}{:<5g}{:<10g}{:<5d}{:<5d}\n'.format(
            self.data['IPOT0'],self.data['FACV0'],self.data['FACW0'],
            self.data['FACVS0'],self.data['FACWS0'],self.data['BETA0'],
            self.data['IMS0'],self.data['IEDG0'])
        out_txt += '{:<5d}{:<5g}{:<5g}{:<5g}{:<5g}{:<10g}{:<5d}{:<5d}\n'.format(
            self.data['IPOT1'],self.data['FACV1'],self.data['FACW1'],
            self.data['FACVS1'],self.data['FACWS1'],self.data['BETA1'],
            self.data['IMS1'],self.data['IEDG1'])
        out_txt += '{:<5d}{:<5g}{:<5g}{:<5g}{:<5g}{:<10g}{:<5d}{:<5d}\n'.format(
            self.data['IPOT2'],self.data['FACV2'],self.data['FACW2'],
            self.data['FACVS2'],self.data['FACWS2'],self.data['BETA2'],
            self.data['IMS2'],self.data['IEDG2'])
        out_txt += '\n\n'+'----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+\n'
        return out_txt 
        
def run_pikoe_from_input_txt(pikoe_input_txt,
                              pikoe_path='../pikoe1.exe',   
                              pikoe_input_path='_test.cnt',
                              verbose=True):
    """
    run pikoe by using input text. 
    Path information must be absolute path!!
    """
    #----need to delete fort before calculation
    import glob 
    from pathlib import Path 
    try: # remove old files 
        os.remove(pikoe_input_path)    
    except:
        if verbose :
            print('Error removing tempoary files')
    print("Be careful that the pikoe_path and input_path to be absolute or cwd")
    print("input path will be cwd for relative path. ")
    #----Assume pikoe file and cnt file are located at the root folder of pikoe 
    root_path = Path(pikoe_input_path).parent 
    
    print('writing input file')
    ff= open(str(pikoe_input_path),'w')
    ff.write(pikoe_input_txt)
    ff.close() 
    print('running pikoe at {}'.format(str(root_path)) )
    #---------
    proc = Popen(pikoe_path +" < "+ str(Path(pikoe_input_path).name)
                 ,stdout=subprocess.PIPE,
                 cwd = str(root_path)
                 ,shell=True)
    (output, err) = proc.communicate()
    proc.wait() 
    print(output);print(err);
    proc.terminate() 
    return output,err           
#    if chck_pikoe_out(fname=pikoe_output_path)==0:
#        out = get_fresco_result()
#        return out 
#    else: 
#        print('There is an Error!!. Check input and output')
#        return 

#=============================================================
if __name__ == '__main__' :
     #--example run Fresco 
     test = PIKOE_input()
     test.set_data(COMMENT='hello') # replace test 
     output = test.write_txt()  
     #------------test using existing sample------------
     #ff = open('pikoe1/sample1/12Cp2pTDXnorm.cnt')
     #output = ff.read()
     #ff.close() 
     #----run 
     out = run_pikoe_from_input_txt(output,
           pikoe_path ="C:/Users/User/Documents/GitHub/GUI/Fresco_GUI/pikoe1/pikoe_window.exe", 
           pikoe_input_path="C:/Users/User/Documents/GitHub/GUI/Fresco_GUI/pikoe1/sample1/_test.cnt",
           verbose=True )    
     
                
