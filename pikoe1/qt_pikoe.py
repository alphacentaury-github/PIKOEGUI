# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 18:16:55 2025

@author: User
"""

import sys
import os
from pathlib import Path 
import webbrowser

from PyQt5.QtWidgets import (QApplication,QMainWindow,
                             QDialog,QFileDialog,QWidget,
                             QTextBrowser,
                             QComboBox,QLabel,QLineEdit,QCheckBox,
                             QMenu,QMenuBar,QDialogButtonBox,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
                             QStackedLayout,
                             QGroupBox,QToolBox,QTabWidget,
                             QPushButton,QTextBrowser,
                             QSpacerItem,QMessageBox,
                             QRadioButton)
from PyQt5 import uic
from PyQt5.QtCore import QSize
from PyQt5.QtCore import Qt
#from PyQt5.QtCore import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

from subprocess import (call,Popen)
import numpy as np

import qt_myutil
from qt_myutil import (combined_Widgets_horizontal,combined_Widgets_vertical,
                       combined_Widgets_grid,QLabel_aligned,
                       text_Browser,WidgetMatplot)
import myutil 
from run_pikoe import *

#==============main window with menu===========================================
class MyWindow(QMainWindow, uic.loadUiType("qt_pikoe_main.ui")[0]):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        size = QSize(900,500)
        self.resize(size)
        self.setWindowTitle('PIKOE GUI')
        # data for main window
        self.gui_setting_path='.pikoe_gui'
        self.path_data = myutil.dict_files(pikoe_path='./pikoe1/pikoe_windows.exe')
        self.load_path_info()
        #---connect Menu action
        self.actionSave.triggered.connect(self.save_file)
        self.actionOpen.triggered.connect(self.open_file)
        self.actionPath_pikoe.triggered.connect(lambda: self.set_path('pikoe_path'))
        self.actionAbout.triggered.connect(self.show_about)
        self.actionDocumentation.triggered.connect(self.show_documents)
        self.actionBug_Report.triggered.connect(self.show_bugreport)
        self.actionOpen_EXFOR.triggered.connect(lambda: webbrowser.open('https://www-nds.iaea.org/exfor/'))  
        #---add Widgets
        self.pikoe = pikoe_GUI() 
        self.pikoe.pikoe_path = self.path_data.data['pikoe_path']
        self.verticalLayout.addWidget(self.pikoe)
        self.show()  
        self.check_pikoe_path() 

    def show_documents(self,):
        webbrowser.open('https://www.sciencedirect.com/science/article/pii/S0010465523004034',new=2 )
        return 

    def show_bugreport(self,):
        webbrowser.open('https://www.rcnp.osaka-u.ac.jp/~kazuyuki/pikoe/index.php',new=2 )
        return 


    def load_path_info(self,):
        # try load path information 
        try:
            self.path_data.load_from_file(self.gui_setting_path)
        except:
            # if file does not exist, create
            self.path_data.save_to_file(self.gui_setting_path)
        
    def check_pikoe_path(self,):
        #---check pikoe_path 
        file_path = Path(self.path_data.data['pikoe_path'])
        if file_path.is_file():
            pass 
        else:
            QMessageBox.warning(self, 'Warning: pikoe path check', 
                f"'{file_path}' does not exist. Please set pikoe location in Configure Menu")            
            print(f"'{file_path}' is not an existing file or does not exist.")
        return 

    def set_path(self,option='pikoe_path'):
        """
        set path for executables 
        possibly later change it to cover other path
        """
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
                      "QFileDialog.getOpenFileName()",
                      "","All Files (*);;Python Files (*.py)",
                      options=options)
        if fileName:
            self.path_data.put_values(**{ option:  fileName})
            self.path_data.save_to_file(self.gui_setting_path) #save gui settings 
            if option=='pikoe_path': #--executable of pikoe path 
                self.pikoe.pikoe_path = fileName 

    def save_file(self,):
        """
        To save/load all reaction calculation
        """
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getSaveFileName(self,
                    "Save file",
                    "",
                    ",All Files (*);",
                    options=options)
        para_dict   = self.pikoe.get_values() 
        if fileName:
            #>> json
            jj = json.dumps(para_dict,indent=2)
            ff = open(fileName,'w')
            ff.write(jj)
            ff.close()
            return para_dict
        else :
            print('No file is chosen')
            return

    def open_file(self,):
        """
        To save/load all reaction calculation
        """
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getOpenFileName(self,
                    "Open file",
                    "",
                    ",All Files (*);",
                    options=options)
        if fileName:
            #>> json
            with open(fileName,"r") as ff:
                jj = json.load(ff)
            self.pikoe.put_values(jj)
            return jj
        else :
            print('No file is chosen')
            return
        return 


    def show_about(self,):
        xxx = QLabel(
            "PIKOE GUI\n"+
            "written by Y.-H.Song(IRIS,IBS), K. Yoshida(RCNP)\n"
            )
        dialog = qt_myutil.CustomDialog(xxx) 
        result = dialog.exec_()
        
        return
    
#==============================================================================
class pikoe_GUI(QWidget,):
    def __init__(self,):
        super().__init__() 
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        
        self.pikoe_input = PIKOE_input() # input texts  
        self.pikoe_path = '' #executable file path 
        
        #----list of files to run pikoe 
        self.pikoe_files = {
            'input'  : [100,'unknown','_test.cnt'], # runtime input # this is the reference path 
            'KIBTBL' : [10,'unknown','tbl_12Cp2pTDXnorm.dat'], # table output
            'KIBOUT' : [6,'unknown','12Cp2pTDXnorm.outlist'], # basic information output 
            'KIBELM' : [11,'old','elem/nnampFL.dat'], #input elementary 
            'ISH': [12,'old','in_ish.dat'],    # ISH>9 case, input  
            'KIBBS': [13,'unknown','out_kibbs.dat'],# KIBBS >0 ,output 
            'IVAR': [14,'old','in_ivar.dat'],  # input, IVAR >9
                    #unit should be in (14,24,34,44)
            'KIBTMD': [15,'unknown','out_kibtmd.dat'], # output KIBTMD > 0       
            'KIBLG': [16,'unknown','out_kiblg.dat'], # output KIBTMD > 0       
            'KIBPX': [17,'unknown','out_kibpx.dat'], # output KIBPX > 0       
            'KIBTR': [18,'unknown','out_kibtr.dat'], # output KIBTR > 0       
            'KIBTL': [19,'unknown','out_kibtl.dat'], # output KIBTL > 0       
            'IPOT0': [20,'old','pot/EDAD1p12C_e.dat'], # input ipot0 >9 
            'IPOT1': [21,'old','pot/EDAD1p11B_e.dat'], # input ipot1 >9 
            'IPOT2': [22,'old','pot/EDAD1p11B_e.dat'] # input ipot2 >9 
            }
        
        #------------------------------------
        # create all input widgets 
        #------------------------------------
        self.widget_dict={} 
        
        self.widget_dict['COMMENT'] = QLineEdit('12C(p,2p)11B_gs@392MeV DWIA TDX normal')
        self.widget_dict['LIMFS'] = QLineEdit('1000')
        self.widget_dict['LIMFS'].setToolTip(
            """limfs       Limit of size of output file (unit # = kibtbl) in unit of MB
              If 0, replaced by 1000000 (roughly 1TB)
              One line is estimated to be 128 bytes."""
            )
        
        self.widget_dict['IONS'] = QComboBox()
        self.widget_dict['IONS'].addItems(['0 : No ','1 : Yes'])
        self.widget_dict['IONS'].setCurrentIndex(0)
        self.widget_dict['IONS'].setToolTip(
            """  ions=0      No output for the configuration that is kinematically not allowed
      <>0     Output for the configuration that is kinematically not allowed"""
            )
        
        self.widget_dict['IFRM'] = QComboBox()
        self.widget_dict['IFRM'].addItems(['0 : Lab','1 : CM','2: Projectile-rest' ])
        self.widget_dict['IFRM'].setCurrentIndex(0) 
        self.widget_dict['IFRM'].setToolTip(
            """  ifrm=0      Kinematics/observable output given in the Laboratory frame (L)
      =1                                            the c.m. frame (G)
      =2                                            the Projectile-rest frame (V)"""
            )
                
        self.widget_dict['IMIR'] = QComboBox()
        self.widget_dict['IMIR'].addItems(['0 : default','1 : reverse z'])
        self.widget_dict['IMIR'].setCurrentIndex(0) 
        self.widget_dict['IMIR'].setToolTip(
            """  imir<>0     +z -> -z conversion in kinematics output""" )
        
        self.widget_dict['ICAL'] = QComboBox()
        self.widget_dict['ICAL'].addItems(['0 : Survey Mode','1 :Computing Mode'])
        self.widget_dict['ICAL'].setCurrentIndex(1) 
        self.widget_dict['ICAL'].setToolTip(
            """  ical=0      Kinematics survey mode; observables not calculated
      <>0     Observables calculated    """
            )
        
        self.widget_dict['ZP'] = QLineEdit('1.0')
        self.widget_dict['AP'] = QLineEdit('1.007825')
        self.widget_dict['ZA'] = QLineEdit('6.0')
        self.widget_dict['AA'] = QLineEdit('12.0')
        self.widget_dict['ZP'].setToolTip(
            """zp ap     Atomic # and mass (in u) of projectile""")
        self.widget_dict['AP'].setToolTip(
            """zp ap     Atomic # and mass (in u) of projectile""")
        self.widget_dict['ZA'].setToolTip(
            """za aa     Atomic # and mass (in u) of target""")
        self.widget_dict['AA'].setToolTip(
            """zp ap     Atomic # and mass (in u) of target""")
        
        self.widget_dict['IKIN'] = QComboBox()
        self.widget_dict['IKIN'].addItems(['0: normal','1: inverse'])
        self.widget_dict['IKIN'].setCurrentIndex(0)
        self.widget_dict['IKIN'].setToolTip(
            " ikin=0      Normal kinematics (probe is the projectile)\n"
           +"    <>0     Inverse kinematics (nucleus to be broken is the projectile)\n"  ) 
        
        self.widget_dict['ELAB'] = QLineEdit('392.0')
        self.widget_dict['ELAB'].setToolTip(
            "elab        Kinetic energy per nucleon of the projectile (in MeV)"
            )
        
        self.widget_dict['ICTREIN'] = QComboBox()
        self.widget_dict['ICTREIN'].addItems(['0 : ','1 : '])
        self.widget_dict['ICTREIN'].setCurrentIndex(0)
        self.widget_dict['ICTREIN'].setToolTip(
          """  ictrein=0   The total kinetic energy Ein of the projectile is calculated with\n
              Ein = elab*nint(ap) or elab*nint(aa).
         <>0  Ein = elab*ap or elab*aa is used."""
            )
        
        self.widget_dict['ISH'] = QComboBox() 
        self.widget_dict['ISH'].addItems(['0: fixed Depth',
                                          '1: Adjust Depth',
                                          '11: External file'])
        self.widget_dict['ISH'].setCurrentIndex(1)
        self.widget_dict['ISH'].setToolTip(
          """  ish=0       Depth (in MeV, >0) of s.p. pot. is specified by ebind.
     =1       Depth of s.p. pot. is changed to reproduce the binding energy
              given by ebind (in MeV, >0).
     >9       S.p. wave function read from external file (unit # = ish)"""  
            )
        
        
        self.widget_dict['EBIND'] = QLineEdit('15.96')
        self.widget_dict['EBIND'].setToolTip(
          """  ebind       Binding energy (in MeV, >0) of the struck nucleon (when ish=1 or ish>9)
              When ish=0 and ibmc=0, regarded as the depth (in MeV) of the central
              s.p. pot.
              When ish=0 and ibmc<>0, it has no meaning."""  
            )
        
        self.widget_dict['ZSP'] = QLineEdit('1.0')
        self.widget_dict['ASP'] = QLineEdit('1.007825')
        self.widget_dict['ZSP'].setToolTip(
            """zsp asp     Atomic # and mass (in u) of particle N to be struck""")
        self.widget_dict['ASP'].setToolTip(
            """zsp asp     Atomic # and mass (in u) of particle N to be struck""")
        
        self.widget_dict['BETASP'] =QLineEdit('0.85')
        self.widget_dict['BETASP'].setToolTip(
            """  betasp      Range of nonlocality (in fm) for particle N in A (Perey-Buck correction)
              If negative, the correction function is read from the external file
              (unit # = ish); if ish<10, process is terminated."""
            )
        
        self.widget_dict['ICTRM'] = QComboBox() 
        self.widget_dict['ICTRM'].addItems(['1: particle mass','2: reduced mass'])
        self.widget_dict['ICTRM'].setCurrentIndex(0)
        self.widget_dict['ICTRM'].setToolTip(
           """  ictrm       Mass parameter in the Sch. Eq. for the bound-state problem:
              =1  Mass of particle N in A
               2  Reduced mass""" 
            )
        
        self.widget_dict['FJ'] = QLineEdit('1.5')
        self.widget_dict['FL'] = QLineEdit('1.0')
        self.widget_dict['FJ'].setToolTip(
           """fj fl       Total and orbital angular momenta of the s.p. state""")
        self.widget_dict['FL'].setToolTip(
           """fj fl       Total and orbital angular momenta of the s.p. state""")
        
        self.widget_dict['SFAC'] = QLineEdit('1.77')
        self.widget_dict['SFAC'].setToolTip(
           """  sfac        Spectroscopic factor (maximum = 2*fj+1 for nucleon-knockout)""")
        self.widget_dict['NOD'] = QLineEdit('0')
        self.widget_dict['NOD'].setToolTip(
            """nod         # of nodes (0 for the lowest state)""")
        
        self.widget_dict['KIBBS'] = QComboBox() # QLineEdit('0') #change into ComboBox !! 
        self.widget_dict['KIBBS'].addItems(['0 : no','1: output'])
        self.widget_dict['KIBBS'].setCurrentIndex(0)
        self.widget_dict['KIBBS'].setToolTip(
            """  kibbs       Unit # for the bound state wave function (BSWF) output
              If kibbs > 0, R, BSWF(R), nonlocality correction function, and
              corrected BSWF (normalized) are written.""")
        
        
        self.widget_dict['IBMC'] = QComboBox() #central pot 
        self.widget_dict['IBMC'].addItems(['1: BM pot','2: Parameter input'])
        self.widget_dict['IBMC'].setCurrentIndex(0)
        self.widget_dict['IBMC'].setToolTip(
            """  ibmc=1      Bohr-Mottelson (BM) pot. used for central part
              If ish>9 and betasp>0, ibmc must be 1.
      <>1     Parameters read-in for central part"""
            )
        
        self.widget_dict['RC'] = QLineEdit('1.35')
        self.widget_dict['RC'].setToolTip(
            """rc          Reduced radial parameter (in fm)""")
        self.widget_dict['RC'].setEnabled(False)
        
        self.widget_dict['ICTRC'] = QComboBox()
        self.widget_dict['ICTRC'].addItems(['0: aa^1/3','1: (aa-asp)^1/3',
                                            '2: (aa-asp)^1/3 + asp^1/3','3: 1'])
        self.widget_dict['ICTRC'].setCurrentIndex(1)
        self.widget_dict['ICTRC'].setToolTip(
            """  ictrc       Control of rc
              =0: fac=aa^1/3
               1: fac=(aa-asp)^1/3
               2: fac=(aa-asp)^1/3 + asp^1/3
               3: fac=1
              The radial parameter is given by rc*fac.""")
        self.widget_dict['ICTRC'].setEnabled(False)      
        
        self.widget_dict['A0C'] = QLineEdit('0.65')
        self.widget_dict['A0C'].setToolTip(
            """a0c         Diffuseness parameter (in fm)""")
        self.widget_dict['A0C'].setEnabled(False)
        
        self.widget_dict['RCL'] = QLineEdit('1.35') 
        self.widget_dict['RCL'].setToolTip(
            """rcl         Reduced Coulomb radius (in fm)"""
            )
        self.widget_dict['RCL'].setEnabled(False) 
        
        self.widget_dict['ICTRCL'] = QComboBox()
        self.widget_dict['ICTRCL'].addItems(['0: aa^1/3','1: (aa-asp)^1/3',
                                             '2: (aa-asp)^1/3 + asp^1/3','3: 1'])
        self.widget_dict['ICTRCL'].setCurrentIndex(1)
        self.widget_dict['ICTRCL'].setToolTip(
            """ictrcl      Control of rcl (similar to ictrc)"""
            )
        self.widget_dict['ICTRCL'].setEnabled(False) 

      

        self.widget_dict['IBMS'] = QComboBox() #spin-orbit pot
        self.widget_dict['IBMS'].addItems(['1: BM pot','2: Parameter input'])
        self.widget_dict['IBMS'].setCurrentIndex(0)
        self.widget_dict['IBMS'].setToolTip(
           """  ibms=1      BM pot. used for spin-orbit part
      <>1     Parameters read-in for spin-orbit part""" 
            )

        self.widget_dict['V0LS'] = QLineEdit('8.2')
        self.widget_dict['V0LS'].setToolTip(
            'v0ls        Depth (in MeV, >0) of spin orbit part')
        self.widget_dict['V0LS'].setEnabled(False)
        
        self.widget_dict['RS'] = QLineEdit('1.35')
        self.widget_dict['RS'].setToolTip(
            """ rs          Reduced radial parameter (in fm)""")
        self.widget_dict['RS'].setEnabled(False)
        
        self.widget_dict['ICTRS'] = QComboBox()
        self.widget_dict['ICTRS'].addItems(['0: aa^1/3','1: (aa-asp)^1/3',
                                            '2: (aa-asp)^1/3 + asp^1/3','3: 1'])
        self.widget_dict['ICTRS'].setCurrentIndex(1) 
        self.widget_dict['ICTRS'].setToolTip(
            """ictrs       Control of rs (similar to ictrc)"""
            )
        self.widget_dict['ICTRS'].setEnabled(False)
        
        self.widget_dict['AS'] = QLineEdit('0.65')
        self.widget_dict['AS'].setToolTip(
            'as          Diffuseness parameter (in fm)'
            )
        self.widget_dict['AS'].setEnabled(False) 
        
        
        
        
        self.widget_dict['LMAX0'] = QLineEdit('-90')
        self.widget_dict['LMAX0'].setToolTip(
            """  lmaxi       Maximum orbital angular momentum for particle i (=0,1,2)
              If negative, lmaxi=min( nint(K_i*R_max), |lmaxi| ) (automatic setting)""")
        self.widget_dict['LMAX0'].setEnabled(False)
        
        self.widget_dict['LMAX1'] = QLineEdit('-90')
        self.widget_dict['LMAX1'].setEnabled(False)
        self.widget_dict['LMAX2'] = QLineEdit('-90')
        self.widget_dict['LMAX2'].setEnabled(False) 
        self.widget_dict['LMAX1'].setToolTip(
            """  lmaxi       Maximum orbital angular momentum for particle i (=0,1,2)
              If negative, lmaxi=min( nint(K_i*R_max), |lmaxi| ) (automatic setting)""")
        self.widget_dict['LMAX2'].setToolTip(
            """  lmaxi       Maximum orbital angular momentum for particle i (=0,1,2)
              If negative, lmaxi=min( nint(K_i*R_max), |lmaxi| ) (automatic setting)""")

        self.widget_dict['IVAR'] = QComboBox() #kinematics profile.. 
        self.widget_dict['IVAR'].addItems(['1: TDX(T1,Ω1,Ω2)',
                                           '2: TDX(KB,ΩB,Ω2)',
                                           '3: QDX(T1,Ω1,T2,ϕ2)',
                                           '9: ...',
                                           '14: external TDX(T1,Ω1,Ω2)',
                                           '24: external TDX(KB,ΩB,Ω2) ',
                                           '34: external QDX(T1,Ω1,T2,ϕ2)',
                                           '44: external...'
                                           ]
            )
        self.widget_dict['IVAR'].setCurrentIndex(0)
        self.widget_dict['IVAR'].setToolTip(
            """  [NOTE: Units of energy (T), angle (theta and phi), and wave number (K) are MeV, degree,
         and fm^-1, respectively. These quantities must be given in the L-frame.]
  ivar=1      T1, theta1, phi1, theta2, and phi2 controlled in
              L11, L12, L13, L14, and L15, respectively
              TDX with respect to T1, Omega1, and Omega2 calculated
      =2      K_B, thetaB, phiB, theta2, and phi2 controlled in
              L11, L12, L13, L14, and L15, respectively
              TDX with respect to K_B, OmegaB, and Omega2 calculated
      =3      T1, theta1, phi1, T2, and phi2 controlled in
              L11, L12, L13, L14, and L15, respectively
              QDX with respect to T1, Omega1, T2, phi2 calculated
      =9      K_Bz in the A-frame controlled in L11 and
              K_Bb in the A-frame controlled in L12
              L13-L15 have no meaning.
              MDs of nucleus B in the A-frame calculated
      >9      Kinematics profile in the L-frame read from external file (unit # = ivar)
              L11-L15 have no meaning.
              (t1,th1,ph1,t2,th2,ph2,kb,thb,phb) to be given with 9f10.0
               9 < ivar < 20: t1,th1,ph1,th2,ph2 used to evaluate the other 4
                              TDX with respect to T1, Omega1, and Omega2 calculated
              19 < ivar < 30: kb,thb,phb,th2,ph2 used to evaluate the other 4
                              TDX with respect to K_B, OmegaB, and Omega2 calculated
                              unit of kb controlled by kunt (see below)
              29 < ivar < 40: t1,th1,ph1,t2,ph2 used to evaluate the other 4
                              QDX with respect to T1, Omega1, T2, phi2 calculated
              39 < ivar     : all of 9 used and energy-momentum cons. checked 
                              TDX with respect to T1, Omega1, and Omega2 calculated
                              unit of kb controlled by kunt (see below)"""
            )
        
        self.widget_dict['IEX'] = QComboBox() 
        self.widget_dict['IEX'].addItems(['0: default','1: ...'])
        self.widget_dict['IEX'].setCurrentIndex(0) 
        self.widget_dict['IEX'].setToolTip(
            """  iex<>1      Particle 1 is particle 0
     =1       Particle 1 is particle N in A in the initial state (not applicable when
              ielm=4 or ielm=6)"""
            )
        
        self.widget_dict['FKNCUT'] = QLineEdit('2.5')
        self.widget_dict['FKNCUT'].setToolTip(
            """fkncut      Cutoff for the missing momentum (in fm^-1)"""
            )
        
        
        self.widget_dict['IXUNT'] = QComboBox() 
        self.widget_dict['IXUNT'].addItems(['1: micro barn','2: milli barn'])
        self.widget_dict['IXUNT'].setCurrentIndex(0) 
        self.widget_dict['IXUNT'].setToolTip(
            """ixunt       Control for the unit of cross section in output
            =1 : micro barn
            <>1: milli barn"""
            )
        
        self.widget_dict['KUNT'] = QComboBox()
        self.widget_dict['KUNT'].addItems(['0: fm^-1','1: MeV/c','2: GeV/c'])
        self.widget_dict['KUNT'].setCurrentIndex(1)
        self.widget_dict['KUNT'].setToolTip( 
            """  kunt        Control for the unit of K_B in output
              Also applied to external file for kinematics profile (unit # = ivar)
              =0: fm^-1
              =1: MeV/c
              =2: GeV/c"""
            )
        
        
        self.widget_dict['IVVAR'] = QComboBox()   #QLineEdit('0')
        self.widget_dict['IVVAR'].addItems(['0: fixed at min','1: ranges'])
        self.widget_dict['IVVAR'].setCurrentIndex(0)
        self.widget_dict['VARMIN'] = QLineEdit('251.0')
        self.widget_dict['VARMAX'] = QLineEdit('255.0')
        self.widget_dict['VARMAX'].setEnabled(False)
        self.widget_dict['DVAR'] = QLineEdit('10.0')
        self.widget_dict['DVAR'].setEnabled(False) 
        for item in ['IVVAR','VARMIN','VARMAX','DVAR']:
            self.widget_dict[item].setToolTip(
            """  ivx=0       x (= var, thx, phx, et2, or ph2) is fixed at xmin
     <>0      x varies from xmin to xmax with the increment dx""")
     

        self.widget_dict['IVTHX'] = QComboBox()   #QLineEdit('0')
        self.widget_dict['IVTHX'].addItems(['0: fixed at min','1: ranges'])
        self.widget_dict['IVTHX'].setCurrentIndex(0) 
        self.widget_dict['THXMIN'] = QLineEdit('32.5')
        self.widget_dict['THXMAX'] = QLineEdit('180.0')
        self.widget_dict['THXMAX'].setEnabled(False)
        self.widget_dict['DTHX'] = QLineEdit('10.0')
        self.widget_dict['DTHX'].setEnabled(False) 
        for item in ['IVTHX','THXMIN','THXMAX','DTHX']:
            self.widget_dict[item].setToolTip(
            """  ivx=0       x (= var, thx, phx, et2, or ph2) is fixed at xmin
     <>0      x varies from xmin to xmax with the increment dx""")
                

        self.widget_dict['IVPHX'] = QComboBox() #QLineEdit('0')
        self.widget_dict['IVPHX'].addItems(['0: fixed at min','1: ranges']  )
        self.widget_dict['IVPHX'].setCurrentIndex(0)
        self.widget_dict['PHXMIN'] = QLineEdit('0.0')
        self.widget_dict['PHXMAX'] = QLineEdit('40.0')
        self.widget_dict['PHXMAX'].setEnabled(False)
        self.widget_dict['DPHX'] = QLineEdit('10.0')
        self.widget_dict['DPHX'].setEnabled(False) 
        for item in ['IVPHX','PHXMIN','PHXMAX','DPHX']:
            self.widget_dict[item].setToolTip(
            """  ivx=0       x (= var, thx, phx, et2, or ph2) is fixed at xmin
     <>0      x varies from xmin to xmax with the increment dx""")     
        

        self.widget_dict['IVET2'] = QComboBox()  #QLineEdit('1')
        self.widget_dict['IVET2'].addItems(['0: fixed at min','1: ranges']  )
        self.widget_dict['IVET2'].setCurrentIndex(1)
        self.widget_dict['ET2MIN'] = QLineEdit('0.0')
        self.widget_dict['ET2MAX'] = QLineEdit('180.0')
        self.widget_dict['DET2'] = QLineEdit('0.5')
        for item in ['IVET2','ET2MIN','ET2MAX','DET2']:
            self.widget_dict[item].setToolTip(
            """  ivx=0       x (= var, thx, phx, et2, or ph2) is fixed at xmin
     <>0      x varies from xmin to xmax with the increment dx""")     
     

        self.widget_dict['IVPH2'] = QComboBox() #QLineEdit('0')
        self.widget_dict['IVPH2'].addItems(['0: fixed at min','1: ranges']  )
        self.widget_dict['IVPH2'].setCurrentIndex(0) 
        self.widget_dict['PH2MIN'] = QLineEdit('180.0')
        self.widget_dict['PH2MAX'] = QLineEdit('360.0')
        self.widget_dict['PH2MAX'].setEnabled(False)
        self.widget_dict['DPH2'] = QLineEdit('10.0')
        self.widget_dict['DPH2'].setEnabled(False) 
        for item in ['IVPH2','PH2MIN','PH2MAX','DPH2']:
            self.widget_dict[item].setToolTip(
            """  ivx=0       x (= var, thx, phx, et2, or ph2) is fixed at xmin
     <>0      x varies from xmin to xmax with the increment dx""")     
     

        #----file specifications
        self.widget_dict['KIBTBL'] = QComboBox()
        self.widget_dict['KIBTBL'].addItems(['10 : output table'])
        self.widget_dict['KIBTBL'].setToolTip(
            """  kibtbl      Unit # for TDX table when ivar<>9
              Unit # for momentum distribution in cylindrical representation with phiB=0
              when ivar=9"""
            )
        
        self.widget_dict['KIBOUT'] = QComboBox()
        self.widget_dict['KIBOUT'].addItems(['6 : outlist'])
        self.widget_dict['KIBOUT'].setToolTip(
            """kibout      Unit # for outlist file"""
            )
        
        self.widget_dict['KIBTMD'] = QComboBox()
        self.widget_dict['KIBTMD'].addItems(['0 : no','15 : write'])
        self.widget_dict['KIBTMD'].setCurrentIndex(0)
        self.widget_dict['KIBTMD'].setToolTip(
            """  kibtmd      Unit # for transition matrix density (when > 0)
              Allowed only when single or no degree of freedom in kinematics is varied.
              Not allowed when ielm=4 or 6.""" 
            )
        self.widget_dict['KIBLG'] = QComboBox() 
        self.widget_dict['KIBLG'].addItems(['0 : no','16 : out'])
        self.widget_dict['KIBLG'].setCurrentIndex(0)
        self.widget_dict['KIBLG'].setToolTip(
            """kiblg       Unit # for longitudinal momentum distribution (when > 0)"""
            )
        self.widget_dict['KIBPX'] = QComboBox() # QLineEdit('0')
        self.widget_dict['KIBPX'].addItems(['0: no','17 : out'])
        self.widget_dict['KIBPX'].setCurrentIndex(0)
        self.widget_dict['KIBPX'].setToolTip(
            """kibpx       Unit # for p_x momentum distribution (when > 0)
              if > 0, ivthx=1 and (thxmax-thxmin)/dthx+1 > 3 must be satisfied"""
            )        
        self.widget_dict['KIBTR'] = QComboBox() #QLineEdit('0')
        self.widget_dict['KIBTR'].addItems(['0 : no','18 : out'])
        self.widget_dict['KIBTR'].setCurrentIndex(0)
        self.widget_dict['KIBTR'].setToolTip(
            """kibtr       Unit # for transverse momentum distribution (when > 0)"""
            )
        
        self.widget_dict['KIBTL'] = QComboBox() #QLineEdit('0')
        self.widget_dict['KIBTL'].addItems(['0 : no','19: out']   )
        self.widget_dict['KIBTL'].setCurrentIndex(0) 
        self.widget_dict['KIBTL'].setToolTip(
            """  kibtl       Unit # for total momentum distribution (when > 0)
              if > 0, ivthx=1, (thxmax-thxmin)/dthx+1 > 3,
                      ivvar=1, and (varmax-varmin)/dvar+1 > 3 must be satisfied
[NOTE: kibtbl and kibout must be positive. kiblg, kibpx, kibtr, kibtl are effective
       only when ivar=9. Trivially, all unit #s must be different from each other.]"""
            )
 
        self.widget_dict['IELM'] = QComboBox() # QLineEdit('3')
        self.widget_dict['IELM'].addItems(['0 : sigma',
                                           '3 : dsigma',
                                           '4 : t-matrix'])
        self.widget_dict['IELM'].setCurrentIndex(2) 
        self.widget_dict['IELM'].setToolTip(
            """  ielm=0      Isotropic free NN cross section at elab (in mb)
      =3      Free differential NN cross section (in mb/sr) read from external file
              (unit # = kibelm) 
      =4      Free (on-shell) NN t-matrix read from external file (unit # = kibelm) 
              Available only when reaction is in coplanar (see below)"""
            )
      
        self.widget_dict['KIBELM'] = QComboBox() #QLineEdit('11')
        self.widget_dict['KIBELM'].addItems(['11 : ...'])
        self.widget_dict['KIBELM'].setCurrentIndex(0)
        self.widget_dict['KIBELM'].setToolTip(
            """kibelm      Unit # of external file for elementary process (needed when ielm=3 or 4)"""
            )
        
        self.widget_dict['IONSH'] = QComboBox() 
        self.widget_dict['IONSH'].addItems(['1: Final-state','2: Initial-state',
                                            '3: Energy-average','4: Momentum-average'])
        self.widget_dict['IONSH'].setCurrentIndex(0)
        self.widget_dict['IONSH'].setToolTip(
            """  ionsh       Choice of on-shell approximation (needed when ielm=3 or 4)     
       =1     Final-state prescription
       =2     Initial-state prescription
       =3     Energy-average prescription
       =4     Momentum-average prescription"""
            )
        
        self.widget_dict['KINELM'] = QComboBox()   #QLineEdit('0')
        self.widget_dict['KINELM'].addItems(['0 : no  output','1 : output']  )
        self.widget_dict['KINELM'].setCurrentIndex(0)
        self.widget_dict['KINELM'].setToolTip(
            """  kinelm=1    Output for the kinematics of the elementary process
              Always set to 0 when ivar=9
        <>1   No output for the kinematics of the elementary process"""
            )
        
        self.widget_dict['IELMEDG'] = QComboBox()   # QLineEdit('1') 
        self.widget_dict['IELMEDG'].addItems(['0 : ','1: ']  ) 
        self.widget_dict['IELMEDG'].setCurrentIndex(1)
        self.widget_dict['IELMEDG'].setToolTip(
            """  ielmedg=0   If the scattering energy is out of the range prepared, 
              the process is terminated.
         <>0  If the scattering energy is out of the range prepared,
              the value at the nearest energy (on an edge) is adopted."""
            )
        #----
        self.widget_dict['RMAX'] = QLineEdit('15.0')
        self.widget_dict['RMAX'].setToolTip(
            """rmax dr     Maximum and increment of radial coordinate R (in fm)""")
        self.widget_dict['DR'] = QLineEdit('0.1')
        self.widget_dict['DR'].setToolTip(
            """rmax dr     Maximum and increment of radial coordinate R (in fm)""")
        self.widget_dict['NGR'] = QLineEdit('30')
        self.widget_dict['NGTH'] = QLineEdit('30')
        self.widget_dict['NGPH'] = QLineEdit('40')
        self.widget_dict['NGK1'] = QLineEdit('0')
        self.widget_dict['NGPH1Q'] = QLineEdit('0')
        for item in ['NGR','NGTH','NGPH','NGK1','NGPH1Q']:
            self.widget_dict[item].setToolTip(
            """     ngx      # of nodes for Gauss-Legendre quadratures for x=R, theta_R, phi_R, K1,
              or phi_1Q
              ng for K1 and phi_1Q effective only when ivar=9"""
            )
        
        self.widget_dict['IPOT0'] = QComboBox()
        self.widget_dict['IPOT0'].addItems(['0 : plane wave',
                                            '1 : KD pot',
                                            '20: external file'])
        self.widget_dict['IPOT0'].setCurrentIndex(2)        
        self.widget_dict['FACV0'] = QLineEdit('1.0')
        self.widget_dict['FACW0'] = QLineEdit('1.0')
        self.widget_dict['FACVS0'] = QLineEdit('1.0')
        self.widget_dict['FACWS0'] = QLineEdit('1.0')
        self.widget_dict['BETA0'] = QLineEdit('-0.85')
        self.widget_dict['IMS0'] = QLineEdit('0')
        self.widget_dict['IEDG0'] = QLineEdit('1')
        
        self.widget_dict['IPOT1'] = QComboBox()
        self.widget_dict['IPOT1'].addItems(['0 : plane wave',
                                            '1 : KD pot',
                                            '21: external file'])
        self.widget_dict['IPOT1'].setCurrentIndex(2)
        self.widget_dict['FACV1'] = QLineEdit('1.0')
        self.widget_dict['FACW1'] = QLineEdit('1.0')
        self.widget_dict['FACVS1'] = QLineEdit('1.0')
        self.widget_dict['FACWS1'] = QLineEdit('1.0')
        self.widget_dict['BETA1'] = QLineEdit('-0.85')
        self.widget_dict['IMS1'] = QLineEdit('0')
        self.widget_dict['IEDG1'] = QLineEdit('1')
        
        self.widget_dict['IPOT2'] = QComboBox()
        self.widget_dict['IPOT2'].addItems(['0 : plane wave',
                                            '1 : KD pot',
                                            '22: external file'])
        self.widget_dict['IPOT2'].setCurrentIndex(2)
        self.widget_dict['FACV2'] = QLineEdit('1.0')
        self.widget_dict['FACW2'] = QLineEdit('1.0')
        self.widget_dict['FACVS2'] = QLineEdit('1.0')
        self.widget_dict['FACWS2'] = QLineEdit('1.0')
        self.widget_dict['BETA2'] = QLineEdit('-0.85')
        self.widget_dict['IMS2'] = QLineEdit('0')
        self.widget_dict['IEDG2'] = QLineEdit('1')
        
        for item in ['0','1','2']:
            self.widget_dict['IPOT'+item].setToolTip( 
        """  ipoti=0     Plane wave is used for particle i (=0,1,2).
        1     Koning-Delaroche pot. (with Coulomb) is adopted for particle i.
       >9     Optical pot. for particle i is read from external file (unit # = ipoti);
              it must cover the energy range appropriately.""")
            self.widget_dict['FACV'+item].setToolTip(
            """      
  facvi       Real part of central term multiplied by facvi
              If negative, the Coulomb potential is turned off
              and |facvi| is used as facvi.""")
            self.widget_dict['FACW'+item].setToolTip(
            """
  facwi       Imaginary part of central term multiplied by facwi""")
            self.widget_dict['FACVS'+item].setToolTip(
            """
  facvsi      Real part of spin-orbit term multiplied by facvsi""")
            self.widget_dict['FACWS'+item].setToolTip(
            """
  facwsi      Imaginary part of spin-orbit term multiplied by facwsi""")
            self.widget_dict['BETA'+item].setToolTip(
            """
  betai       Range of nonlocality (in fm) for particle i.
              If negative, the correction function is read from the external file
              (unit # = ipoti); if ipoti<10, process is terminated.
            """)
            self.widget_dict['IMS'+item].setToolTip(
            """
  imsi=0      Reduced energy used as the kinematical factor
      <>0     Reduced mass used as the kinematical factor """)
            self.widget_dict['IEDG'+item].setToolTip(""" 
  iedgi=0     If the scattering energy is out of the range prepared,
              the process is terminated.
       <>0    If the scattering energy is out of the range prepared,
              the value at the nearest energy (on an edge) is adopted.""")
        
        #-----actions-----------------------------------------------------------
 
        self.widget_dict['IVAR'].currentIndexChanged.connect( 
            self.ivar_changed )
        self.widget_dict['IVVAR'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IVVAR'].currentIndex()==0, 
                [self.widget_dict['VARMAX'],self.widget_dict['DVAR']])
            ) 
        self.widget_dict['IVTHX'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IVTHX'].currentIndex()==0, 
                [self.widget_dict['THXMAX'],self.widget_dict['DTHX']])
            ) 
        
        self.widget_dict['IVPHX'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IVPHX'].currentIndex()==0, 
                [self.widget_dict['PHXMAX'],self.widget_dict['DPHX']])
            ) 
        self.widget_dict['IVET2'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IVET2'].currentIndex()==0, 
                [self.widget_dict['ET2MAX'],self.widget_dict['DET2']])
            ) 
        self.widget_dict['IVPH2'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IVPH2'].currentIndex()==0, 
                [self.widget_dict['PH2MAX'],self.widget_dict['DPH2']])
            ) 
        self.widget_dict['ISH'].currentIndexChanged.connect( self.ish_changed )
        self.widget_dict['KIBBS'].currentIndexChanged.connect(lambda: self.set_filename('KIBBS'))      
        self.widget_dict['IBMC'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IBMC'].currentIndex()==0, 
                [self.widget_dict['RC'],self.widget_dict['ICTRC'],
                 self.widget_dict['A0C'],self.widget_dict['RCL'],
                 self.widget_dict['ICTRCL'] ])
            )   
        self.widget_dict['IBMS'].currentIndexChanged.connect( 
            lambda : self.change_enables(self.widget_dict['IBMS'].currentIndex()==0, 
                [self.widget_dict['V0LS'],self.widget_dict['RS'],
                 self.widget_dict['ICTRS'],self.widget_dict['AS']])
            ) 
        self.widget_dict['IELM'].currentIndexChanged.connect(self.ielm_changed)      
        self.widget_dict['IPOT0'].currentIndexChanged.connect(lambda: self.ipot_changed('IPOT0'))
        self.widget_dict['IPOT1'].currentIndexChanged.connect(lambda: self.ipot_changed('IPOT1'))
        self.widget_dict['IPOT2'].currentIndexChanged.connect(lambda: self.ipot_changed('IPOT2'))
        #---------Buttons------------------------------------------------------
        self.dict_buttons = {}
        
        self.dict_buttons['RUN'] = QPushButton('Run Pikoe')
        self.dict_buttons['RUN'].clicked.connect(self.pikoe_input_run)
        self.dict_buttons['RUN'].setStyleSheet("background-color: green;")
        
        self.dict_buttons['input'] = QPushButton('set cnt file')
        self.dict_buttons['input'].clicked.connect(
            lambda : self.set_filename(target_var_name='input') ) 
        self.dict_buttons['input'].setToolTip('reference path in pikoe run')
        
        self.dict_buttons['KIBTBL'] = QPushButton('set KIBTBL file')
        self.dict_buttons['KIBTBL'].clicked.connect(
            lambda : self.set_filename(target_var_name='KIBTBL') ) 
        
        self.dict_buttons['KIBOUT'] = QPushButton('set KIBOUT file')
        self.dict_buttons['KIBOUT'].clicked.connect(
            lambda : self.set_filename(target_var_name='KIBOUT') ) 
                
        self.dict_buttons['READ_KIBTBL'] = QPushButton('Read KIBTBL')
        self.dict_buttons['READ_KIBTBL'].clicked.connect(
            lambda: self.load_pikoe_output('KIBTBL') )
        self.dict_buttons['READ_KIBOUT'] = QPushButton('READ_KIBOUT')
        self.dict_buttons['READ_KIBOUT'].clicked.connect(
            lambda: self.load_pikoe_output('KIBOUT') )
        
        self.dict_buttons['PLOT_KIBTBL'] = QPushButton('Plot KIBTBL')
        self.dict_buttons['PLOT_KIBTBL'].clicked.connect(
            #lambda: QMessageBox.warning(self,'Warning','Not available yet'  ) 
            self.plot_widget_matplot
            )
        #----------------------------------------------------------------------
        #-------display widgets------------------------------------------------
        #----------------------------------------------------------------------
                                   
        #-----------combined widgets-------------------------
        self.widg_L1 = combined_Widgets_horizontal(
                         [QLabel('Comment:'),self.widget_dict['COMMENT']])
        self.widg_L2 = combined_Widgets_horizontal(
                         [QLabel('LIMFS:'),self.widget_dict['LIMFS'],
                          QLabel('IONS'),self.widget_dict['IONS'],
                          QLabel('IFRM'),self.widget_dict['IFRM'],
                          QLabel('IMIR'),self.widget_dict['IMIR'],
                          QLabel('ICAL'),self.widget_dict['ICAL'] ])
        self.widg_L3L4 = combined_Widgets_horizontal([QLabel('Projectile: Charge'),self.widget_dict['ZP'],
                                     QLabel('mass(amu)'),self.widget_dict['AP'],
                                     QLabel('Taget: charge'),self.widget_dict['ZA'],
                                     QLabel('mass(amu)'),self.widget_dict['AA'],
                                     QLabel('kinematics:'),self.widget_dict['IKIN'],
                                     QLabel('elab(MeV/u)'),self.widget_dict['ELAB']
                                     ])
        self.widg_L5L6 = combined_Widgets_horizontal([
            QLabel('K.O. particle: charge'),self.widget_dict['ZSP'],
            QLabel(' mass(amu)'),self.widget_dict['ASP'],
            QLabel('s.p. wave option: '),self.widget_dict['ISH'],
            QLabel('B.E.(MeV)'),self.widget_dict['EBIND'],
            QLabel('w.f. : node'),self.widget_dict['NOD'],
            QLabel('L'),self.widget_dict['FL'],
            QLabel('J'),self.widget_dict['FJ'],
            QLabel('kibbs'),self.widget_dict['KIBBS']
            ])
        self.widg_L7 = combined_Widgets_horizontal(
            [QLabel('Central pot. :'),self.widget_dict['IBMC'],
             QLabel('reduced radius: r'), self.widget_dict['RC'],
             QLabel('factor'), self.widget_dict['ICTRC'],
             QLabel('diffuseness'), self.widget_dict['A0C'],
             QLabel('reduced Coulomb: r'), self.widget_dict['RCL'],
             QLabel('factor'), self.widget_dict['ICTRCL']])
        
        self.widg_L8 = combined_Widgets_horizontal(
            [QLabel('LS pot. :'),self.widget_dict['IBMS'],
             QLabel('Depth(>0)'),self.widget_dict['V0LS'],
             QLabel('reduced radius: r'), self.widget_dict['RS'],
             QLabel('factor'), self.widget_dict['ICTRS'],
             QLabel('diffuseness'), self.widget_dict['AS']])
        
        self.L_checkbox = QCheckBox('Default') 
        self.L_checkbox.toggle() 
        self.L_checkbox.stateChanged.connect(self.change_L_check)
        
        self.widg_L9 = combined_Widgets_horizontal(
            [QLabel('Maximum L:'), self.L_checkbox,  
             QLabel('L0'), self.widget_dict['LMAX0'],
             QLabel('L1'), self.widget_dict['LMAX1'],
             QLabel('L2'), self.widget_dict['LMAX2'],
             ]
            )
        self.widg_L10 = combined_Widgets_horizontal( 
            [QLabel('Kinematics profile:'),self.widget_dict['IVAR'],  
             QLabel('iex'),self.widget_dict['IEX'],
             QLabel('Cutoff(fm^-1)'),self.widget_dict['FKNCUT'],
             QLabel('xs. unit'),self.widget_dict['IXUNT'],
             QLabel('mom. unit'),self.widget_dict['KUNT']]
            )
        
        
        self.widg_L11_15 = combined_Widgets_grid( 
           [ [QLabel('ivvar:'),self.widget_dict['IVVAR'],  
             QLabel('min'),self.widget_dict['VARMIN'],
             QLabel('max'),self.widget_dict['VARMAX'],
             QLabel('dx'),self.widget_dict['DVAR']],
            [QLabel('ivthx:'),self.widget_dict['IVTHX'],  
             QLabel('min'),self.widget_dict['THXMIN'],
             QLabel('max'),self.widget_dict['THXMAX'],
             QLabel('dx'),self.widget_dict['DTHX'] ],
            [QLabel('ivphx:'),self.widget_dict['IVPHX'],  
             QLabel('min'),self.widget_dict['PHXMIN'],
             QLabel('max'),self.widget_dict['PHXMAX'],
             QLabel('dx'),self.widget_dict['DPHX']],
            [QLabel('ivet2:'),self.widget_dict['IVET2'],  
              QLabel('min'),self.widget_dict['ET2MIN'],
              QLabel('max'),self.widget_dict['ET2MAX'],
              QLabel('dx'),self.widget_dict['DET2']],
            [QLabel('ivph2:'),self.widget_dict['IVPH2'],  
             QLabel('min'),self.widget_dict['PH2MIN'],
             QLabel('max'),self.widget_dict['PH2MAX'],
             QLabel('dx'),self.widget_dict['DPH2'] ] 
           ],opt='json' )
        

        self.widg_L16 = combined_Widgets_horizontal(
            [QLabel('kib: tbl'),self.widget_dict['KIBTBL'],
             QLabel('out'),self.widget_dict['KIBOUT'],
             QLabel('tmd'),self.widget_dict['KIBTMD'],
             QLabel('lg'),self.widget_dict['KIBLG'],
             QLabel('px'),self.widget_dict['KIBPX'],
             QLabel('tr'),self.widget_dict['KIBTR'],
             QLabel('tl'),self.widget_dict['KIBTL']]
            )         
        self.widg_L17 = combined_Widgets_horizontal(
            [QLabel('ielm/kibelm'),self.widget_dict['IELM'],
             #QLabel('kibelm'),self.widget_dict['KIBELM'],
             QLabel('ionsh'),self.widget_dict['IONSH'],
             QLabel('kinelm'),self.widget_dict['KINELM'],
             QLabel('ielmedg'),self.widget_dict['IELMEDG']]
            )
        self.widg_L18 = combined_Widgets_horizontal(
            [QLabel('Rmax'),self.widget_dict['RMAX'],
             QLabel('dr'),self.widget_dict['DR'],
             QLabel('ngr'),self.widget_dict['NGR'],
             QLabel('ngth'),self.widget_dict['NGTH'],
             QLabel('ngph'),self.widget_dict['NGPH'],
             QLabel('ngk1'),self.widget_dict['NGK1'],
             QLabel('ngph1q'),self.widget_dict['NGPH1Q']] )
        self.widg_L19_21 = combined_Widgets_grid(
          [  [QLabel('ipot0'),self.widget_dict['IPOT0'],
             QLabel('facv0'),self.widget_dict['FACV0'],
             QLabel('facw0'),self.widget_dict['FACW0'],
             QLabel('facvs0'),self.widget_dict['FACVS0'],
             QLabel('facws0'),self.widget_dict['FACWS0'],
             QLabel('beta0'),self.widget_dict['BETA0'],
             QLabel('ims0'),self.widget_dict['IMS0'],
             QLabel('iedg0'),self.widget_dict['IEDG0']],
             [QLabel('ipot1'),self.widget_dict['IPOT1'],
             QLabel('facv1'),self.widget_dict['FACV1'],
             QLabel('facw1'),self.widget_dict['FACW1'],
             QLabel('facvs1'),self.widget_dict['FACVS1'],
             QLabel('facws1'),self.widget_dict['FACWS1'],
             QLabel('beta1'),self.widget_dict['BETA1'],
             QLabel('ims1'),self.widget_dict['IMS1'],
             QLabel('iedg1'),self.widget_dict['IEDG1']],
             [QLabel('ipot2'),self.widget_dict['IPOT2'],
              QLabel('facv2'),self.widget_dict['FACV2'],
              QLabel('facw2'),self.widget_dict['FACW2'],
              QLabel('facvs2'),self.widget_dict['FACVS2'],
              QLabel('facws2'),self.widget_dict['FACWS2'],
              QLabel('beta2'),self.widget_dict['BETA2'],
              QLabel('ims2'),self.widget_dict['IMS2'],
              QLabel('iedg2'),self.widget_dict['IEDG2']] 
          ],opt='json')
        
        self.widg_files = combined_Widgets_horizontal(
            [self.dict_buttons['input'], 
             self.dict_buttons['KIBTBL'] ,
             self.dict_buttons['KIBOUT']
             ])
     
        #-----put into tabs---------------------------------- 
        self.tabs = QTabWidget()
        
        self.tab1 = combined_Widgets_vertical(
            [self.widg_L1, self.widg_L2,self.widg_L3L4,self.widg_L16,
             self.widg_files]
            )
        self.tab2 = combined_Widgets_vertical(
            [self.widg_L5L6 , self.widg_L7,self.widg_L8]
            )
        self.tab3 = combined_Widgets_vertical(
            [self.widg_L10,self.widg_L11_15]
            )
        self.tab4 = combined_Widgets_vertical(
            [self.widg_L9,self.widg_L17,self.widg_L18,
             self.widg_L19_21]
            )
        self.tabs.addTab(self.tab1,'Basic Setup' )
        self.tabs.addTab(self.tab2,'Structure' )
        self.tabs.addTab(self.tab3,'Kinematics' )
        self.tabs.addTab(self.tab4,'Reaction' )

                                 
        
        self.widg_after_run = combined_Widgets_horizontal(
            [self.dict_buttons['READ_KIBOUT'],
             self.dict_buttons['READ_KIBTBL'],
             self.dict_buttons['PLOT_KIBTBL']
             ])           
        
        #------Final display-------------------------------
        self.layout.addWidget(QLabel(
            'PIKOE: Proton Induced KnockOut reaction calculation for Exclusive process'))
        self.layout.addWidget(self.tabs)
        self.layout.addWidget(self.dict_buttons['RUN'] )
        self.layout.addWidget(self.widg_after_run)
                      
    #==============Methods===============================================    
    def check_pikoe_filepath(self,file_list):
        #---check files exists  file_list is in pikoe_header form 
        check = 0 
        for ff in file_list:
           file_path = Path(ff[2])
           if file_path.is_file():
               pass  
           else:
               if ff[1] in ['old']:  
                   QMessageBox.warning(self, 'Warning: file path check', 
               f"'{file_path}' does not exist. Check file location.")
               check = 1     
        return check 
    
    def change_L_check(self,state):
        if state == Qt.Checked:
            self.widget_dict['LMAX0'].setEnabled(False)
            self.widget_dict['LMAX0'].setText('-90')
            self.widget_dict['LMAX1'].setEnabled(False)
            self.widget_dict['LMAX1'].setText('-90')
            self.widget_dict['LMAX2'].setEnabled(False)
            self.widget_dict['LMAX2'].setText('-90')
        else :
            self.widget_dict['LMAX0'].setEnabled(True)
            self.widget_dict['LMAX1'].setEnabled(True)
            self.widget_dict['LMAX2'].setEnabled(True)
        return
    
    def change_enables(self,condition,widgets_list):
        """ 
        switch enables status of list of widgets by state 
        state==0 : disabled 
             ==1 : enabled 
        """
        if condition :
            for item in widgets_list:
                item.setEnabled(False)
        else:
            for item in widgets_list:
                item.setEnabled(True)
        return 
        
    def ivar_changed(self,):
        string =  self.widget_dict['IVAR'].currentText() 
        ivar_val = int(string.split(':')[0])
        widgets_L11_12 = [ self.widget_dict['IVVAR'],  
          self.widget_dict['VARMIN'],
          self.widget_dict['VARMAX'],
          self.widget_dict['DVAR'],
          self.widget_dict['IVTHX'],  
          self.widget_dict['THXMIN'],
          self.widget_dict['THXMAX'],
          self.widget_dict['DTHX'] ] 
        widgets_L13_15 = [  self.widget_dict['IVPHX'],  
          self.widget_dict['PHXMIN'],
          self.widget_dict['PHXMAX'],
          self.widget_dict['DPHX'],
          self.widget_dict['IVET2'],  
          self.widget_dict['ET2MIN'],
          self.widget_dict['ET2MAX'],
          self.widget_dict['DET2'],
          self.widget_dict['IVPH2'],  
          self.widget_dict['PH2MIN'],
          self.widget_dict['PH2MAX'],
          self.widget_dict['DPH2'] ] 
        if ivar_val in [1,2,3]:
            for item in widgets_L11_12:
                item.setEnabled(True) 
            for item in widgets_L13_15:
                item.setEnabled(True)  
               
        elif ivar_val==9:
            for item in widgets_L11_12:
                item.setEnabled(True) 
            for item in widgets_L13_15:
                item.setEnabled(False)  
                
        elif ivar_val in [14,24,34,44]: #external file 
            for item in widgets_L11_12:
                item.setEnabled(False) 
            for item in widgets_L13_15:
                item.setEnabled(False)  
            self.set_filename(target_var_name='IVAR')
        return 
 
    
     
    
    
 
    def ish_changed(self,):
        string =  self.widget_dict['ISH'].currentText() 
        ish_val = int(string.split(':')[0])
        if ish_val >9 :
            self.set_filename(target_var_name='ISH')            
        return 

    def ielm_changed(self,):
        #--note that this is related with KIBELM 
        string =  self.widget_dict['IELM'].currentText() 
        val = int(string.split(':')[0])
        if val==0 :
            pass
        elif val in [3,4]:    
            self.set_filename(target_var_name='KIBELM')            
        return 

    def ipot_changed(self,target_var_name='IPOT0'):
        string =  self.widget_dict[target_var_name].currentText() 
        val = int(string.split(':')[0])
        if val > 9 : 
            self.set_filename(target_var_name)
        return     

    def set_filename(self,target_var_name=''):
        """ 
        get file path for target_var from user 
        need to be relative path to pikoe executable file 
        """
        options = QFileDialog.Options()
           
        fstate = self.pikoe_files[target_var_name][1] 
        default_filename = self.pikoe_files[target_var_name][2] 
        default_filename = str(Path(self.pikoe_path).parent) +'/'+default_filename
        #----get new file path 
        if fstate in ['new','unknown']:
            fileName, _filter = QFileDialog.getSaveFileName(self,
                    "set filename",
                    default_filename,
                    "All Files (*)",
                    options=options)
            #Cancel returns None 
            #if file exist, ask whether to change it. 
            #  --> if cancel, at this moment, no change occurs. 
        #----open old file 
        elif fstate in ['old']:
            fileName, _filter = QFileDialog.getOpenFileName(self,
                    "set filename",
                    default_filename,
                    "All Files (*)",
                    options=options)
            print(fileName)
        #---store the data in pikoe
        #   because of length limit, it have to be relative path 
        if fileName:
            try:
                with open(fileName,'x') as f:
                    f.write('')
            except FileExistsError:
                pass 
            
            #method1
            # base = Path(self.pikoe_path).parent  
            # target = Path(fileName) 
            # rel_path = target.relative_to(base) if base in target.parents or target == base else target 
            # rel_path = f"{rel_path.as_posix()}" #fortran format         
            # self.pikoe_files[target_var_name][2] = './'+rel_path       
            #method2
            base_dir = os.path.dirname(self.pikoe_path)
            target_dir = os.path.dirname(fileName)
            target_name = os.path.basename(fileName)
            rel_path = Path(os.path.relpath(target_dir,start=base_dir))
            rel_path = f"{rel_path.as_posix()}" #fortran format
            
            self.pikoe_files[target_var_name][2] = rel_path+'/'+target_name  
            print(self.pikoe_files[target_var_name][2])
            
        return self.pikoe_files[target_var_name][2]
        
            
    def get_values(self,):
        self.data = {}         
        for key in self.widget_dict.keys():
            widg = self.widget_dict[key]
            index = key 
            if isinstance(widg, QLabel):
                self.data[index] = widg.text()
            elif isinstance(widg, QLineEdit):
                self.data[index] = widg.text()
            elif isinstance(widg, QComboBox):
                #----get actual pikoe input value 
                self.data[index] = widg.currentText().split(':')[0]
            elif isinstance(widg, QRadioButton):
                self.data[index] = widg.isChecked()
            elif isinstance(widg, QCheckBox):
                self.data[index] = widg.isChecked()     
            else: #other widgets are skipped
                pass 
            
        return self.data

    def put_values(self,data):
        # inverse of get_values
        # key is number with get_values. 0, 1
        # however, json makes the key as string, '0','1'
        for key in self.widget_dict.keys():
            widg = self.widget_dict[key]
            index = key 
            if isinstance(widg, QLabel):
                widg.setText(data[index]) # all data are strings
            elif isinstance(widg, QLineEdit):
                widg.setText(data[index])
            elif isinstance(widg, QComboBox):
                #--convert pikoe input value to index 
                str_ = data[index] 
                for i in range(widg.count() ):
                    if str_ in widg.itemText(i).split(':')[0]:
                        break  
                widg.setCurrentIndex(i)
            elif isinstance(widg,QRadioButton):
                widg.setChecked(data[index])
            elif isinstance(widg, QCheckBox):
                widg.setChecked(data[index])
            else: #other widgets are skipped
                pass
            
    def get_files(self,):
        """ 
        convert GUI pikoe_files into file control form of pikoe input 
        """
        file_list = []
        for key in self.pikoe_files.keys(): #file related items 
            if (key in ['ISH','IVAR','IPOT0','IPOT1','IPOT2']
                and int(self.data[key]) > 9 ): #(value>9)
                self.pikoe_files[key][0] = int(self.data[key]) #unit numbers from gui
                file_list.append(self.pikoe_files[key])
            elif ( key in ['KIBTBL','KIBOUT','KIBELM',
                          'KIBTMD','KIBPX','KIBTR','KIBTL','KIBBS'] 
                  and int(self.data[key]) != 0 ):
                self.pikoe_files[key][0] = int(self.data[key]) #unit numbers from gui
                file_list.append(self.pikoe_files[key])
            elif (key in ['input','output']):
                pass  
        return file_list         
            
            
    def pikoe_input_run(self,):
        """ 
        convert GUI inputs input pikoe input text file 
        and run the pikoe. 
        """
        self.get_values() 
        print('Warning! check the index/actual input consistency!!')
        self.pikoe_input.set_data(**self.data) 
        file_list = self.get_files() # file data settings     
        check = self.check_pikoe_filepath(file_list)
        self.pikoe_input.set_data(HEADERS = file_list)
        if check!=0: 
            QMessageBox.about(self,'pikoe run','Check file paths')
            return         
        self.check_pikoe_filepath(file_list)
        
        txt_input = self.pikoe_input.write_txt() 
        print(txt_input) 
        QMessageBox.about(self,'pikoe run','Starting Pikoe. Please press Ok and wait')
        output,err = run_pikoe_from_input_txt(txt_input,pikoe_path=self.pikoe_path,
               pikoe_input_path=self.pikoe_files['input'][2])
        QMessageBox.about(self,'pikoe run',
                output.decode("utf-8") +'\n'
                +'Calculation is done\n')
        return 
    
    def load_pikoe_output(self,target_var):
        filename = self.pikoe_files[target_var][2] 
        with open(filename,'r') as ff:
            txt = ff.read() 
        xxx = QDialog() 
        vbox = QVBoxLayout() 
        xxx.setLayout(vbox)
        tb  = QTextBrowser()
        vbox.addWidget(tb)
        tb.append(txt)
        size = QSize(600,500)
        xxx.resize(size)
        xxx.exec_() 
        return 
    
    def plot_widget_matplot(self,):
        """ 
        Plot KIBTBL testing...
        incomplete .... 
        """
        filename = self.pikoe_files['KIBTBL'][2] 
        data = np.loadtxt(filename,skiprows=1)
        
        mm = QWidget() 
        mm_layout = QVBoxLayout()
        mm.setLayout(mm_layout)
        
        widg_plot = qt_myutil.WidgetMatplot()
        mm_layout.addWidget(widg_plot)
        #---------------------------
        #-----test plot KIBTBL data 
        #---------------------------
        btn1 = QPushButton('Plot')
        btn2 = QPushButton('Exp. data')
        
        colx = QComboBox()
        colx.addItems(['0 : t1l[MeV]','1: th1l[deg]',
                       '2 : ph1l[deg]','3 : t2l[MeV]',
                       '4 : th2l[deg]','5 : ph2l[deg]',
                       '6 : pbl[MeV/c]','7: thbl[deg]',
                       '8 : phbl[deg]','9: pr[MeV/c]',
                       '10: isol','11: tdx[ub/(MeVsr2)]','12 : Ay'])
        colx.setCurrentIndex(4)
        coly = QComboBox()
        coly.addItems(['0 : t1l[MeV]','1: th1l[deg]',
                       '2 : ph1l[deg]','3 : t2l[MeV]',
                       '4 : th2l[deg]','5 : ph2l[deg]',
                       '6 : pbl[MeV/c]','7: thbl[deg]',
                       '8 : phbl[deg]','9: pr[MeV/c]',
                       '10: isol','11: tdx[ub/(MeVsr2)]','12 : Ay'])
        coly.setCurrentIndex(12)
        
        widg_xy = combined_Widgets_horizontal([QLabel('Column x'),colx,
                                    QLabel('Column y'),coly,btn2,btn1])
        
        
        mm_layout.addWidget(widg_xy)
        
        btn1.clicked.connect(lambda: btn1_clicked() )
        btn2.clicked.connect(lambda: QMessageBox.about(
                mm,'exp. data',
                'Not available function yet.') )
        
        def btn1_clicked():
            fig = plt.Figure()
            ax = fig.add_subplot(111)
            x = colx.currentIndex()
            y = coly.currentIndex()
            ax.plot(data[:,x],data[:,y])
            widg_plot.rm_plot() 
            widg_plot.add_plot(fig) 
            return     
        
        dialog = qt_myutil.CustomDialog(mm)
        dialog.exec_() 
        return 
    
#===================================================================================================
if __name__ == "__main__":
    def run_app():
        """
        launcher of Qt in Spyder
        to avoid error
        """
        if not QApplication.instance():
            app = QApplication(sys.argv)
        else:
            app = QApplication.instance()

        #myWindow = pikoe_GUI() #without menu 
        myWindow = MyWindow() # with menu 
        myWindow.show()
        app.exec_()
        return myWindow
    m = run_app()        