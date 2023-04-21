#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Print ORFeus logo
# Author: Mary Richardson
# Date: 2023.03.31
# -----------------------------------------------------------------------------


CYAN = '\033[96m'
GREY = '\033[97m'
END = '\033[0m'

def logo():

    print('\n'+ \
          '                                                                     '+ CYAN+'//\\\\'+END+'\n'+ \
          '                                                                     '+ CYAN+'\\\\//'+END+'\n'+ \
          '                                                                  '+ CYAN+'//\\\\'+END+'\n'+ \
          '                                                                  '+ CYAN+'\\\\//'+END+'\n'+ \
          '                                                                     '+ CYAN+'//\\\\'+END+'\n'+ \
          '                                                                     '+ CYAN+'\\\\//'+END+'\n'+ \
          '                                                                        '+ CYAN+'//\\\\'+END+'\n'+ \
          '                                                                        '+ CYAN+'\\\\//'+END+'\n'+ \
          '                                                                      '+ CYAN+'//\\\\'+END+'\n'+ \
          '                                                                      '+ CYAN+'\\\\//'+END+'\n'+ \
          '                                                                         '+ CYAN+'//\\\\'+END+'\n'+ \
          CYAN+'    ACGACG    GUCGUGC   CAGUCGGCA'+END+ '                                        '+ CYAN+'\\\\//'+END+'\n'+ \
          CYAN+'  AGUCAGUCGU  ACGUUACG  UGACGUAAC'+END+ '                                   '+ GREY+'.RIBOSOMERIBOSOME.'+END+'\n'+ \
          CYAN+' GUC      AUG CA    GUA AGC'+END+ '        CGUAG   CAG   CGA  CAGUGC     '+ GREY+'.RIBOSOMERIBOSOMERIBOSO.'+END+'\n'+ \
          CYAN+' CA        GC UGAGUUGC  GUACUGAA'+END+ ' GUA   AGC ACG   GAC CGU       '+ GREY+'.MERIBOSOMERIBOSOMERIBOSOME.'+END+'\n'+ \
          CYAN+' ACG      CGU GG ACGU   CGAGUAG'+END+ '  ACGUAGUAG UCA   AUG  AGCA    '+ GREY+'.RIBOSOMERIBOSOMERIBOSOMERIBO.'+END+'\n'+ \
          CYAN+'  CGCGGGCAUG  AC   ACG  UGA'+END+ '      GGA       CAU   GCU     GCA  '+ GREY+'.SOMERIBOSOMERIBOSOMERIBOSOME.'+END+'\n'+ \
          CYAN+'    GUCAGC    GC    GCC CAG'+END+ '        AGUAG    CAGUGCA  CAUGCA     '+ GREY+'.RIBOSOMERIBOSOMERIBOSOME.'+END+'\n'+ \
          '                                                                 '+ GREY+'.OMERIBOSOMERIBOSOMERIB.'+END+'\n'+ \
          '##########################################################################################################'+'\n'+ \
          '                                                                 '+ GREY+'.RIBOSOMERIBOSOMERIBOSO.   '+END+'\n'+ \
          '                                                                '+ GREY+'.MERIBOSOMERIBOSOMERIBOSO.  '+END+'\n'+ \
          '                                                                   '+ GREY+'.MERIBOSOMERIBOSOME.     '+END)
