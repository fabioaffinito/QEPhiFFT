	  �  <   k820309              15.0        OU                                                                                                           
       electrons_base.f90 ELECTRONS_BASE                      @                              
       DP                                                                                                                 @@                                                     @                                                      @                                                      @                                                         p          p            p                                    @                                                      @@                                                        p          p            p                                    @                                	                         p          p            p                                    @                                
                      @                                                      @                                                      @                                                         p          p            p                                    @                                                         p          p            p                                    @                                                      @                                                      @@                                                     @                                                         p          p            p                                    @                                                      @                                                     @@                                                 
                &                                                                                          
                @@                                                                  &                                                    @@                                                 
                &                                                    @@                                                                  &                                                    @@                                                                  &                                           #         @                                                   	   #ELECTRONS_BASE_INITVAL%INT    #ELECTRONS_BASE_INITVAL%DBLE    #ELECTRONS_BASE_INITVAL%SUM    #ELECTRONS_BASE_INITVAL%TRIM    #ELECTRONS_BASE_INITVAL%MOD     #ELECTRONS_BASE_INITVAL%MAX !   #ELECTRONS_BASE_INITVAL%NINT "   #ELECTRONS_BASE_INITVAL%ABS #   #ZV_ $   #NA_ %   #NSP_ &   #NBND_ '   #NSPIN_ (   #OCCUPATIONS_ )   #F_INP *   #TOT_CHARGE_ +   #TOT_MAGNETIZATION_ ,                                                                                   INT                                                DBLE                                                SUM                                                TRIM                                                 MOD                                           !     MAX                                           "     NINT                                           #     ABS           
                                 $                   
              &                                                     
                                  %                                 &                                                     
                                  &                     
                                  '                     
                                  (                     
  @                              )                    1           
                                 *                   
              &                   &                                                     
                                 +     
                
  @                              ,     
      #         @                                  -                   #SET_NELUP_NELDW%INT .   #SET_NELUP_NELDW%MOD /   #SET_NELUP_NELDW%NINT 0   #SET_NELUP_NELDW%ABS 1   #TOT_MAGNETIZATION_ 2   #NELEC_ 3   #NELUP_ 4   #NELDW_ 5                                                                                               .     INT                                           /     MOD                                           0     NINT                                           1     ABS           
  @                              2     
                
  @                              3     
                D                                4     
                 D                                5     
       #         @                                   6                    #DEALLOCATE_ELCT%ALLOCATED 7                                             7     ALLOCATED #         @                                   8                   #DISTRIBUTE_BANDS%MOD 9   #NBGRP :   #MY_BGRP_ID ;                                             9     MOD           
@ @                               :                     
@ @                               ;              �   *      fn#fn    �   C   J  KINDS      p       DP+KINDS    }  @       NBND    �  @       NBNDX    �  @       NSPIN    =  �       NEL    �  @       NELT      �       NUPDWN    �  �       IUPDWN    9  @       NUDX    y  @       NBSP    �  @       NBSPX    �  �       NUPDWN_BGRP    �  �       IUPDWN_BGRP    !  @       NUDX_BGRP    a  @       NBSP_BGRP    �  @       NBSPX_BGRP    �  �       I2GUPDWN_BGRP (   u  @       TELECTRONS_BASE_INITVAL    �  @       KEEP_OCC    �  �       F    �  @       QBAC    �  �       ISPIN    M	  �       F_BGRP    �	  �       ISPIN_BGRP    e
  �       IBGRP_G2L '   �
  �      ELECTRONS_BASE_INITVAL +   �  <      ELECTRONS_BASE_INITVAL%INT ,     =      ELECTRONS_BASE_INITVAL%DBLE +   O  <      ELECTRONS_BASE_INITVAL%SUM ,   �  =      ELECTRONS_BASE_INITVAL%TRIM +   �  <      ELECTRONS_BASE_INITVAL%MOD +     <      ELECTRONS_BASE_INITVAL%MAX ,   @  =      ELECTRONS_BASE_INITVAL%NINT +   }  <      ELECTRONS_BASE_INITVAL%ABS +   �  �   a   ELECTRONS_BASE_INITVAL%ZV_ +   E  �   a   ELECTRONS_BASE_INITVAL%NA_ ,   �  @   a   ELECTRONS_BASE_INITVAL%NSP_ -     @   a   ELECTRONS_BASE_INITVAL%NBND_ .   Q  @   a   ELECTRONS_BASE_INITVAL%NSPIN_ 4   �  L   a   ELECTRONS_BASE_INITVAL%OCCUPATIONS_ -   �  �   a   ELECTRONS_BASE_INITVAL%F_INP 3   �  @   a   ELECTRONS_BASE_INITVAL%TOT_CHARGE_ :   �  @   a   ELECTRONS_BASE_INITVAL%TOT_MAGNETIZATION_             SET_NELUP_NELDW $     <      SET_NELUP_NELDW%INT $   X  <      SET_NELUP_NELDW%MOD %   �  =      SET_NELUP_NELDW%NINT $   �  <      SET_NELUP_NELDW%ABS 3     @   a   SET_NELUP_NELDW%TOT_MAGNETIZATION_ '   M  @   a   SET_NELUP_NELDW%NELEC_ '   �  @   a   SET_NELUP_NELDW%NELUP_ '   �  @   a   SET_NELUP_NELDW%NELDW_       g       DEALLOCATE_ELCT *   t  B      DEALLOCATE_ELCT%ALLOCATED !   �  }       DISTRIBUTE_BANDS %   3  <      DISTRIBUTE_BANDS%MOD '   o  @   a   DISTRIBUTE_BANDS%NBGRP ,   �  @   a   DISTRIBUTE_BANDS%MY_BGRP_ID 