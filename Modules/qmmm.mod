	  �1  �   k820309              15.0        �OU                                                                                                           
       qmmm.f90 QMMM              QMMM_CONFIG QMMM_INITIALIZATION QMMM_SHUTDOWN QMMM_UPDATE_POSITIONS QMMM_UPDATE_FORCES QMMM_ADD_MM_FIELD                                                    
                            @                              
       IONODE IONODE_ID STDOUT                      @                              
       WORLD_COMM                      @                              
       MP_BCAST MP_BARRIER MP_ABORT                      @                              
       DP                �                                      u #MP_BCAST_I1    #MP_BCAST_R1 
   #MP_BCAST_C1    #MP_BCAST_Z    #MP_BCAST_ZV    #MP_BCAST_IV !   #MP_BCAST_RV &   #MP_BCAST_CV +   #MP_BCAST_L 0   #MP_BCAST_RM 4   #MP_BCAST_CM 9   #MP_BCAST_IM >   #MP_BCAST_IT C   #MP_BCAST_RT H   #MP_BCAST_LV M   #MP_BCAST_LM R   #MP_BCAST_R4D W   #MP_BCAST_R5D \   #MP_BCAST_CT a   #MP_BCAST_C4D f   #MP_BCAST_C5D k   #         @     @                                               #MSG    #SOURCE    #GID 	                                                                                                                           
                                  	           #         @     @                            
                    #MSG    #SOURCE    #GID                                                   
                 
                                                       
                                             #         @     @                                                #MSG    #SOURCE    #GID                                                                    
                                                       
                                             #         @     @                                               #MP_BCAST_Z%LEN    #MP_BCAST_Z%ICHAR    #MP_BCAST_Z%CHAR    #MSG    #SOURCE    #GID                  @                                LEN               @                                ICHAR               @                                CHAR                                                                1           
                                                       
                                             #         @     @                                               #MP_BCAST_ZV%SIZE    #MP_BCAST_ZV%LEN    #MP_BCAST_ZV%ICHAR    #MP_BCAST_ZV%CHAR    #MSG    #SOURCE    #GID                   @                                SIZE               @                                LEN               @                                ICHAR               @                                CHAR ,                                                                            &                                           1           
                                                       
                                              #         @     @                            !                   #MP_BCAST_IV%SIZE "   #MSG #   #SOURCE $   #GID %                 @                           "     SIZE                                           #                                  &                                                     
                                  $                     
                                  %           #         @     @                            &                   #MP_BCAST_RV%SIZE '   #MSG (   #SOURCE )   #GID *                 @                           '     SIZE                                          (                   
               &                                                     
                                  )                     
                                  *           #         @     @                            +                   #MP_BCAST_CV%SIZE ,   #MSG -   #SOURCE .   #GID /                 @                           ,     SIZE                                          -                                  &                                                     
                                  .                     
                                  /           #         @     @                            0                    #MSG 1   #SOURCE 2   #GID 3                                              1                      
                                  2                     
                                  3           #         @     @                           4                   #MP_BCAST_RM%SIZE 5   #MSG 6   #SOURCE 7   #GID 8                 @                           5     SIZE                                          6                   
 	              &                   &                                                     
                                  7                     
                                  8           #         @     @                            9                   #MP_BCAST_CM%SIZE :   #MSG ;   #SOURCE <   #GID =                 @                           :     SIZE                                          ;                                  &                   &                                                     
                                  <                     
                                  =           #         @     @                            >                   #MP_BCAST_IM%SIZE ?   #MSG @   #SOURCE A   #GID B                 @                           ?     SIZE                                           @                                  &                   &                                                     
                                  A                     
                                  B           #         @     @                            C                   #MP_BCAST_IT%SIZE D   #MSG E   #SOURCE F   #GID G                 @                           D     SIZE                                           E                                  &                   &                   &                                                     
                                  F                     
                                  G           #         @     @                            H                   #MP_BCAST_RT%SIZE I   #MSG J   #SOURCE K   #GID L                 @                           I     SIZE                                          J                   
 
              &                   &                   &                                                     
                                  K                     
                                  L           #         @     @                            M                   #MP_BCAST_LV%SIZE N   #MSG O   #SOURCE P   #GID Q                 @                           N     SIZE                                           O                                  &                                                     
                                  P                     
                                  Q           #         @     @                            R                   #MP_BCAST_LM%SIZE S   #MSG T   #SOURCE U   #GID V                 @                           S     SIZE                                           T                                  &                   &                                                     
                                  U                     
                                  V           #         @     @                            W                   #MP_BCAST_R4D%SIZE X   #MSG Y   #SOURCE Z   #GID [                 @                           X     SIZE                                          Y                   
               &                   &                   &                   &                                                     
                                  Z                     
                                  [           #         @     @                            \                   #MP_BCAST_R5D%SIZE ]   #MSG ^   #SOURCE _   #GID `                 @                           ]     SIZE                                          ^                   
               &                   &                   &                   &                   &                                                     
                                  _                     
                                  `           #         @     @                            a                   #MP_BCAST_CT%SIZE b   #MSG c   #SOURCE d   #GID e                 @                           b     SIZE                                          c                                  &                   &                   &                                                     
                                  d                     
                                  e           #         @     @                            f                   #MP_BCAST_C4D%SIZE g   #MSG h   #SOURCE i   #GID j                 @                           g     SIZE                                          h                                  &                   &                   &                   &                                                     
                                  i                     
                                  j           #         @     @                            k                   #MP_BCAST_C5D%SIZE l   #MSG m   #SOURCE n   #GID o                 @                           l     SIZE                                          m                                  &                   &                   &                   &                   &                                                     
                                  n                     
                                  o                                                      p                      @@                               q                                                       r                       @                               s            #         @                                   t                    #GID u             
                                  u           #         @                                  v                    #ERRORCODE w   #GID x             
                                  w                     
                                  x                                                        y                                                         #         @                                   z                   #QMMM_CONFIG%PRESENT {   #MODE |   #COMM }   #VERBOSE ~   #STEP                  @                           {     PRESENT           
 @                               |                     
 @                               }                     
 @                               ~                     
 @                                          #         @                                   �                    #QMMM_INITIALIZATION%TRIM �                 @                           �     TRIM #         @                                   �                     #         @                                   �                                                                 #         @                                   �                    #FORCE �             
                                 �                   
              &                   &                                           #         @                                   �                        �         fn#fn    �   y   b   uapp(QMMM !   /  @   J  PARALLEL_INCLUDE    o  X   J  IO_GLOBAL    �  K   J  MP_WORLD      ]   J  MP    o  C   J  KINDS     �  �      gen@MP_BCAST+MP    Y  f      MP_BCAST_I1+MP #   �  @   a   MP_BCAST_I1%MSG+MP &   �  @   a   MP_BCAST_I1%SOURCE+MP #   ?  @   a   MP_BCAST_I1%GID+MP      f      MP_BCAST_R1+MP #   �  @   a   MP_BCAST_R1%MSG+MP &   %  @   a   MP_BCAST_R1%SOURCE+MP #   e  @   a   MP_BCAST_R1%GID+MP    �  f      MP_BCAST_C1+MP #     @   a   MP_BCAST_C1%MSG+MP &   K  @   a   MP_BCAST_C1%SOURCE+MP #   �  @   a   MP_BCAST_C1%GID+MP    �  �      MP_BCAST_Z+MP &   p  <      MP_BCAST_Z%LEN+MP=LEN *   �  >      MP_BCAST_Z%ICHAR+MP=ICHAR (   �  =      MP_BCAST_Z%CHAR+MP=CHAR "   '	  L   a   MP_BCAST_Z%MSG+MP %   s	  @   a   MP_BCAST_Z%SOURCE+MP "   �	  @   a   MP_BCAST_Z%GID+MP    �	  �      MP_BCAST_ZV+MP )   �
  =      MP_BCAST_ZV%SIZE+MP=SIZE '   �
  <      MP_BCAST_ZV%LEN+MP=LEN +   *  >      MP_BCAST_ZV%ICHAR+MP=ICHAR )   h  =      MP_BCAST_ZV%CHAR+MP=CHAR #   �  �   a   MP_BCAST_ZV%MSG+MP &   5  @   a   MP_BCAST_ZV%SOURCE+MP #   u  @   a   MP_BCAST_ZV%GID+MP    �  |      MP_BCAST_IV+MP )   1  =      MP_BCAST_IV%SIZE+MP=SIZE #   n  �   a   MP_BCAST_IV%MSG+MP &   �  @   a   MP_BCAST_IV%SOURCE+MP #   :  @   a   MP_BCAST_IV%GID+MP    z  |      MP_BCAST_RV+MP )   �  =      MP_BCAST_RV%SIZE+MP=SIZE #   3  �   a   MP_BCAST_RV%MSG+MP &   �  @   a   MP_BCAST_RV%SOURCE+MP #   �  @   a   MP_BCAST_RV%GID+MP    ?  |      MP_BCAST_CV+MP )   �  =      MP_BCAST_CV%SIZE+MP=SIZE #   �  �   a   MP_BCAST_CV%MSG+MP &   �  @   a   MP_BCAST_CV%SOURCE+MP #   �  @   a   MP_BCAST_CV%GID+MP      f      MP_BCAST_L+MP "   j  @   a   MP_BCAST_L%MSG+MP %   �  @   a   MP_BCAST_L%SOURCE+MP "   �  @   a   MP_BCAST_L%GID+MP    *  |      MP_BCAST_RM+MP )   �  =      MP_BCAST_RM%SIZE+MP=SIZE #   �  �   a   MP_BCAST_RM%MSG+MP &   �  @   a   MP_BCAST_RM%SOURCE+MP #   �  @   a   MP_BCAST_RM%GID+MP      |      MP_BCAST_CM+MP )   �  =      MP_BCAST_CM%SIZE+MP=SIZE #   �  �   a   MP_BCAST_CM%MSG+MP &   d  @   a   MP_BCAST_CM%SOURCE+MP #   �  @   a   MP_BCAST_CM%GID+MP    �  |      MP_BCAST_IM+MP )   `  =      MP_BCAST_IM%SIZE+MP=SIZE #   �  �   a   MP_BCAST_IM%MSG+MP &   A  @   a   MP_BCAST_IM%SOURCE+MP #   �  @   a   MP_BCAST_IM%GID+MP    �  |      MP_BCAST_IT+MP )   =  =      MP_BCAST_IT%SIZE+MP=SIZE #   z  �   a   MP_BCAST_IT%MSG+MP &   6  @   a   MP_BCAST_IT%SOURCE+MP #   v  @   a   MP_BCAST_IT%GID+MP    �  |      MP_BCAST_RT+MP )   2  =      MP_BCAST_RT%SIZE+MP=SIZE #   o  �   a   MP_BCAST_RT%MSG+MP &   +  @   a   MP_BCAST_RT%SOURCE+MP #   k  @   a   MP_BCAST_RT%GID+MP    �  |      MP_BCAST_LV+MP )   '  =      MP_BCAST_LV%SIZE+MP=SIZE #   d  �   a   MP_BCAST_LV%MSG+MP &   �  @   a   MP_BCAST_LV%SOURCE+MP #   0  @   a   MP_BCAST_LV%GID+MP    p  |      MP_BCAST_LM+MP )   �  =      MP_BCAST_LM%SIZE+MP=SIZE #   )  �   a   MP_BCAST_LM%MSG+MP &   �  @   a   MP_BCAST_LM%SOURCE+MP #      @   a   MP_BCAST_LM%GID+MP     M   }      MP_BCAST_R4D+MP *   �   =      MP_BCAST_R4D%SIZE+MP=SIZE $   !  �   a   MP_BCAST_R4D%MSG+MP '   �!  @   a   MP_BCAST_R4D%SOURCE+MP $   "  @   a   MP_BCAST_R4D%GID+MP     ["  }      MP_BCAST_R5D+MP *   �"  =      MP_BCAST_R5D%SIZE+MP=SIZE $   #  �   a   MP_BCAST_R5D%MSG+MP '   $  @   a   MP_BCAST_R5D%SOURCE+MP $   A$  @   a   MP_BCAST_R5D%GID+MP    �$  |      MP_BCAST_CT+MP )   �$  =      MP_BCAST_CT%SIZE+MP=SIZE #   :%  �   a   MP_BCAST_CT%MSG+MP &   �%  @   a   MP_BCAST_CT%SOURCE+MP #   6&  @   a   MP_BCAST_CT%GID+MP     v&  }      MP_BCAST_C4D+MP *   �&  =      MP_BCAST_C4D%SIZE+MP=SIZE $   0'  �   a   MP_BCAST_C4D%MSG+MP '   (  @   a   MP_BCAST_C4D%SOURCE+MP $   D(  @   a   MP_BCAST_C4D%GID+MP     �(  }      MP_BCAST_C5D+MP *   )  =      MP_BCAST_C5D%SIZE+MP=SIZE $   >)  �   a   MP_BCAST_C5D%MSG+MP '   **  @   a   MP_BCAST_C5D%SOURCE+MP $   j*  @   a   MP_BCAST_C5D%GID+MP !   �*  @       IONODE+IO_GLOBAL $   �*  @       IONODE_ID+IO_GLOBAL !   *+  @       STDOUT+IO_GLOBAL $   j+  @       WORLD_COMM+MP_WORLD    �+  Q       MP_BARRIER+MP "   �+  @   a   MP_BARRIER%GID+MP    ;,  `       MP_ABORT+MP &   �,  @   a   MP_ABORT%ERRORCODE+MP     �,  @   a   MP_ABORT%GID+MP    -  p       DP+KINDS    �-  �       QMMM_CONFIG $   .  @      QMMM_CONFIG%PRESENT !   W.  @   a   QMMM_CONFIG%MODE !   �.  @   a   QMMM_CONFIG%COMM $   �.  @   a   QMMM_CONFIG%VERBOSE !   /  @   a   QMMM_CONFIG%STEP $   W/  f       QMMM_INITIALIZATION )   �/  =      QMMM_INITIALIZATION%TRIM    �/  H       QMMM_SHUTDOWN &   B0  t       QMMM_UPDATE_POSITIONS #   �0  S       QMMM_UPDATE_FORCES )   	1  �   a   QMMM_UPDATE_FORCES%FORCE "   �1  H       QMMM_ADD_MM_FIELD 