	    D   k820309              15.0        OU                                                                                                           
       radial_grids.f90 RADIAL_GRIDS              NDMX RADIAL_GRID_TYPE DO_MESH CHECK_MESH HARTREE SERIES WRITE_GRID_ON_FILE READ_GRID_FROM_FILE ALLOCATE_RADIAL_GRID NULLIFY_RADIAL_GRID RADIAL_GRID_COPY gen@DEALLOCATE_RADIAL_GRID                      @                              
       DP                                                       u #DEALLOCATE_RADIAL_GRID_S    #DEALLOCATE_RADIAL_GRID_V    #         @     @X                                               #DEALLOCATE_RADIAL_GRID_S%ASSOCIATED    #GRID                  @                                ASSOCIATED           
D @                                                   #RADIAL_GRID_TYPE    #         @     @X                                                #DEALLOCATE_RADIAL_GRID_V%SIZE    #DEALLOCATE_RADIAL_GRID_V%ASSOCIATED    #GRID 	                 @                                SIZE               @                                ASSOCIATED           
D@                               	                                   &                                           #RADIAL_GRID_TYPE                                                 
                                                                                                                                             �              3500                  @                                '                    #MESH    #R    #R2    #RAB    #SQR    #RM1    #RM2    #RM3    #XMIN    #RMAX    #ZMESH    #DX                 � $                                                             �$                                                          
            &                                                       �$                                         P                 
            &                                                       �$                                         �                 
            &                                                       �$                                         �                 
            &                                                       �$                                         (                
            &                                                       �$                                         p                
            &                                                       �$                                         �                
            &                                                        � $                                         	   
                � $                                        
   
                � $                                           
                � $                                           
   #         @                                                      #DO_MESH%SQRT    #DO_MESH%EXP    #DO_MESH%DBLE    #DO_MESH%LOG    #RMAX    #ZMESH    #XMIN    #DX     #IBOUND !   #GRID "                                 @                                SQRT               @                                EXP               @                                DBLE               @                                LOG           
                                      
                
                                      
                
D                                     
                 
                                       
                
                                  !                     D @                               "                    #RADIAL_GRID_TYPE    #         @                                   #                   #CHECK_MESH%ABS $   #CHECK_MESH%SQRT %   #GRID &                                                         @                           $     ABS               @                           %     SQRT           
                                  &                   #RADIAL_GRID_TYPE    #         @                                   '                   #HARTREE%DBLE (   #K )   #NST *   #MESH +   #GRID ,   #F -   #VH .                                 @                           (     DBLE           
  @                               )                     
                                  *                     
                                  +                     
@ @                               ,                   #RADIAL_GRID_TYPE             
                                 -                    
 	   p          5 � p        r +       5 � p        r +                              D @                              .                    
 
    p          5 � p        r +       5 � p        r +                     #         @                                  /                    #F 0   #R 1   #R2 2   #B 3                                                            0                   
     p          p            p                                                                    1                   
     p          p            p                                                                    2                   
     p          p            p                                    D     �                           3                   
     p           & p         p            p                          #         @                                   4                    #IUNIT 5   #GRID 6             
                                  5                     
                                  6                   #RADIAL_GRID_TYPE    #         @                                   7                    #IUNIT 8   #GRID 9             
                                  8                     D                                 9                    #RADIAL_GRID_TYPE    #         @                                  :                    #GRID ;   #MESH <             
D                                 ;                    #RADIAL_GRID_TYPE              
                                  <           #         @                                  =                    #GRID >             
D                                 >                    #RADIAL_GRID_TYPE    #         @                                   ?                    #X @   #Y A             
                                  @                   #RADIAL_GRID_TYPE              
D @                               A                    #RADIAL_GRID_TYPE       �   &      fn#fn "   �   �   b   uapp(RADIAL_GRIDS    �  C   J  KINDS +   �  |       gen@DEALLOCATE_RADIAL_GRID )   I  {      DEALLOCATE_RADIAL_GRID_S 4   �  C      DEALLOCATE_RADIAL_GRID_S%ASSOCIATED .     ^   a   DEALLOCATE_RADIAL_GRID_S%GRID )   e  �      DEALLOCATE_RADIAL_GRID_V .     =      DEALLOCATE_RADIAL_GRID_V%SIZE 4   @  C      DEALLOCATE_RADIAL_GRID_V%ASSOCIATED .   �  �   a   DEALLOCATE_RADIAL_GRID_V%GRID    %  p       DP+KINDS    �  t       NDMX !   	  �       RADIAL_GRID_TYPE &   �  H   a   RADIAL_GRID_TYPE%MESH #     �   a   RADIAL_GRID_TYPE%R $   �  �   a   RADIAL_GRID_TYPE%R2 %   6  �   a   RADIAL_GRID_TYPE%RAB %   �  �   a   RADIAL_GRID_TYPE%SQR %   ^	  �   a   RADIAL_GRID_TYPE%RM1 %   �	  �   a   RADIAL_GRID_TYPE%RM2 %   �
  �   a   RADIAL_GRID_TYPE%RM3 &     H   a   RADIAL_GRID_TYPE%XMIN &   b  H   a   RADIAL_GRID_TYPE%RMAX '   �  H   a   RADIAL_GRID_TYPE%ZMESH $   �  H   a   RADIAL_GRID_TYPE%DX    :  �       DO_MESH      =      DO_MESH%SQRT    R  <      DO_MESH%EXP    �  =      DO_MESH%DBLE    �  <      DO_MESH%LOG      @   a   DO_MESH%RMAX    G  @   a   DO_MESH%ZMESH    �  @   a   DO_MESH%XMIN    �  @   a   DO_MESH%DX      @   a   DO_MESH%IBOUND    G  ^   a   DO_MESH%GRID    �  �       CHECK_MESH    H  <      CHECK_MESH%ABS     �  =      CHECK_MESH%SQRT     �  ^   a   CHECK_MESH%GRID      �       HARTREE    �  =      HARTREE%DBLE    �  @   a   HARTREE%K    9  @   a   HARTREE%NST    y  @   a   HARTREE%MESH    �  ^   a   HARTREE%GRID      �   a   HARTREE%F    �  �   a   HARTREE%VH      t       SERIES    �  �   a   SERIES%F    �  �   a   SERIES%R      �   a   SERIES%R2    �  �   a   SERIES%B #   S  ]       WRITE_GRID_ON_FILE )   �  @   a   WRITE_GRID_ON_FILE%IUNIT (   �  ^   a   WRITE_GRID_ON_FILE%GRID $   N  ]       READ_GRID_FROM_FILE *   �  @   a   READ_GRID_FROM_FILE%IUNIT )   �  ^   a   READ_GRID_FROM_FILE%GRID %   I  \       ALLOCATE_RADIAL_GRID *   �  ^   a   ALLOCATE_RADIAL_GRID%GRID *     @   a   ALLOCATE_RADIAL_GRID%MESH $   C  R       NULLIFY_RADIAL_GRID )   �  ^   a   NULLIFY_RADIAL_GRID%GRID !   �  V       RADIAL_GRID_COPY #   I  ^   a   RADIAL_GRID_COPY%X #   �  ^   a   RADIAL_GRID_COPY%Y 