!INIT PARAMS:
IG = 1 ;  DE = 0
MD = 36;  DTR = 3.141593 / 180.
ML = 16;  DTR2 = 3.141593 / 180.
MT = 0 ;  G65 = 65.
MTS = 0;

!INPUT SECTION
! - open files ( *.INP, *.OUT, *.SUM)
! - read INP to find max values
! - allocate and init arrays
WRITE(*,*) 'READING INPUT DATA FROM FILE.'
! - read INP to get data
! - for each stack: detect if stack is on top a roof. Set LFLAT and GEPIN == 1 if does.
WRITE(*,*) 'END OF READING INPUT DATA FROM FILE.'
!END OF INPUT SECTION

!GEP VALUES SECTION
WRITE(*,*) 'CALCULATING GEP VALUES.'
! - calc min dist between structures (DO 80)

! GEP STACK HEIGHT CALCULATIONS DETERMINE IF A STACK IS WITHIN A GEP 5L AREA OF INFLUENCE AS STAND ALONE STRUCTURES AND TIERS FOR EVERY QUARTER OF A DEGREE
WRITE(*,*) '  Calculating single tier GEP values.'
! - for every 0.25 deg: (DO 100)
!      + project coordinates of stack and tiers corners
!      + calc tier projected wid,len,hgt,L. 
!      + check if stack is inside SIZ-L5
!      + if stack is on SIZ L5 then call "GPC" (it calculates GEP values)

! GEP STACK HEIGHT CALCULATIONS DETERMINE IF A STACK IS WITHIN A GEP 5L AREA OF INFLUENCE AS COMBINED STRUCTURES AND TIERS FOR EVERY QUARTER OF A DEGREE IDENTIFY TIER GROUPS - EXAMINE FOR COMBINING USE ACTUAL HEIGHTS - EACH GROUP FORMED AROUND FIRST TIER EVERY TIER IS USED AS FIRST OR 'FOCAL' TIER IN SUCCESSION
WRITE(*,*) '  Looking for and calculating any group of tiers GEP values'
WRITE(*,*) '    for a wind flow starting at 0.25 degrees.'
! - for every 0.25 deg: (DO 110)
!      + project coordinates of tiers corners (WHY?) 
!      + calc width, height of "Focal Tier"
!      + check if focal tier could be combined with other:
!          - dist(focal, tier) < max(L_focal,L_tier) => combinable!
!          - save combinable tier in TLIST(C1,:), and TNUM(C1)+=1
! (L:1149) FOR SUFFICIENTLY CLOSE STRUCTURES COMBINE IDENTIFIED STRUCTURES BY GROUPS 
! - for each focal "i" building, and "j" tier (DO 120)
!      + C1=(I-1) * MXTRS + J (calc tier absolute ID)
!      + check if TNUM > 1. J tier absolute ID 
!      + for each combinable tier "T1" (DO 122)  => Creates focal subgroups based on common height
!         - get Heights of focal (HTC) and T1 (HTA). only proceed if HTA < HTC or if T1 is the focal tier. Lets HTA be the common
!         height of the group
!         - initialize focal tier values (xmin,xmax,..date (HTA)
!         - for each combinable tier "T2" (DO 123)
!             - if height of T2 > HTA proceed.
!             - if distance T1-T2 < maxL TL1-TL2: merge! (update: xmin,xmax,ymin,ymax)
!             - then calculate wid,len, and L of merged tiers.

     

!END OF GEP VALUES SECTION




!BUILDING DOWNWASH VALUES SECTION 
!Entiendo que calcula los valores para la tabla de los Stacks
!Repite casi lo mismo que antes
 CALCULATING BUILDING DOWNWASH INPUT VALUES.
  
   Calculating single tier downwash values.
   Calculating group of tiers downwash values.


!END OF BUILDING DOWNWASH VALUES SECTION


