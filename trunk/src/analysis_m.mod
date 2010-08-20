V24 analysis_m
14 analysis_m.F90 S582 0
04/26/2010  19:42:59
use params_m private
use params_m private
enduse
D 33 24 644 456 639 7
D 261 24 1010 192 1009 7
D 267 21 9 2 15 52 0 0 0 0 0
 0 25 3 3 25 25
 0 25 25 3 25 25
S 582 24 0 0 0 8 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 27 0 0 0 0 0 0 analysis_m
S 592 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 595 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 632 3 0 0 0 6 1 1 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 633 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 634 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 635 3 0 0 0 16 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 636 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 637 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 639 25 1 params_m simparameters
R 644 5 6 params_m title1 simparameters
R 645 5 7 params_m title2 simparameters
R 646 5 8 params_m ivol simparameters
R 647 5 9 params_m iquen simparameters
R 648 5 10 params_m iprint simparameters
R 649 5 11 params_m nnbrs simparameters
R 650 5 12 params_m iverlet simparameters
R 651 5 13 params_m nlcx simparameters
R 652 5 14 params_m nlcy simparameters
R 653 5 15 params_m nlcz simparameters
R 654 5 16 params_m nm simparameters
R 655 5 17 params_m nspec simparameters
R 656 5 18 params_m dsp simparameters
R 657 5 19 params_m rpad simparameters
R 658 5 20 params_m rcut simparameters
R 659 5 21 params_m rnear simparameters
R 660 5 22 params_m deltat simparameters
R 661 5 23 params_m temprq simparameters
R 662 5 24 params_m tempsp simparameters
R 663 5 25 params_m zlayer simparameters
R 664 5 26 params_m rqke simparameters
R 665 5 27 params_m press simparameters
R 666 5 28 params_m nloops simparameters
R 667 5 29 params_m nsteps simparameters
R 668 5 30 params_m nprint simparameters
R 669 5 31 params_m ntcm simparameters
R 670 5 32 params_m nchkpt simparameters
R 671 5 33 params_m restart simparameters
R 672 5 34 params_m nose simparameters
R 673 5 35 params_m alternate_quench_md simparameters
R 674 5 36 params_m nout simparameters
R 675 5 37 params_m ntape simparameters
R 676 5 38 params_m dumpx1 simparameters
R 677 5 39 params_m tjob simparameters
R 678 5 40 params_m tfinalise simparameters
R 679 5 41 params_m prevsteps simparameters
R 680 5 42 params_m currentstep simparameters
R 681 5 43 params_m laststep simparameters
R 682 5 44 params_m lastprint simparameters
R 683 5 45 params_m lastchkpt simparameters
R 684 5 46 params_m ntc simparameters
R 685 5 47 params_m strx simparameters
R 686 5 48 params_m uselookup simparameters
R 687 5 49 params_m boxtem simparameters
R 688 5 50 params_m bdel2 simparameters
R 689 5 51 params_m bmass simparameters
S 1003 27 0 0 0 8 1034 582 7513 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 update_thermodynamic_sums
S 1004 27 0 0 0 8 1028 582 7539 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 set_thermodynamic_sums
S 1005 27 0 0 0 8 1031 582 7562 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 get_thermodynamic_sums
S 1006 27 0 0 0 8 1038 582 7585 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 runavs
S 1007 27 0 0 0 8 1044 582 7592 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 rdf
S 1008 27 0 0 0 8 1041 582 7596 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 auto
S 1009 25 0 0 0 261 1 582 7601 c 800000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 thermodynamic_sums
S 1010 5 0 0 0 9 1011 582 7620 800004 0 0 0 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1 1010 0 582 0 0 0 0 spe
S 1011 5 0 0 0 9 1012 582 7635 800004 0 0 8 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1010 1011 0 582 0 0 0 0 ske
S 1012 5 0 0 0 9 1013 582 7650 800004 0 0 16 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1011 1012 0 582 0 0 0 0 ste
S 1013 5 0 0 0 9 1014 582 7665 800004 0 0 24 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1012 1013 0 582 0 0 0 0 sh
S 1014 5 0 0 0 9 1015 582 7679 800004 0 0 32 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1013 1014 0 582 0 0 0 0 sth
S 1015 5 0 0 0 9 1016 582 7694 800004 0 0 40 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1014 1015 0 582 0 0 0 0 sf2
S 1016 5 0 0 0 9 1017 582 7709 800004 0 0 48 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1015 1016 0 582 0 0 0 0 svol
S 1017 5 0 0 0 9 1018 582 7725 800004 0 0 56 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1016 1017 0 582 0 0 0 0 spesq
S 1018 5 0 0 0 9 1019 582 7742 800004 0 0 64 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1017 1018 0 582 0 0 0 0 skesq
S 1019 5 0 0 0 9 1020 582 7759 800004 0 0 72 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1018 1019 0 582 0 0 0 0 stesq
S 1020 5 0 0 0 9 1021 582 7776 800004 0 0 80 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1019 1020 0 582 0 0 0 0 shsq
S 1021 5 0 0 0 9 1022 582 7792 800004 0 0 88 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1020 1021 0 582 0 0 0 0 sthsq
S 1022 5 0 0 0 9 1023 582 7809 800004 0 0 96 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1021 1022 0 582 0 0 0 0 sf2sq
S 1023 5 0 0 0 9 1024 582 7826 800004 0 0 104 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1022 1023 0 582 0 0 0 0 svolsq
S 1024 5 0 0 0 9 1025 582 7844 800004 0 0 112 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1023 1024 0 582 0 0 0 0 svolpe
S 1025 5 0 0 0 267 1 582 7862 800004 0 0 120 0 0 261 0 0 0 0 0 0 0 0 0 0 0 1024 1025 0 582 0 0 0 0 sb0
S 1026 6 4 0 0 261 1 582 7877 80003c 0 0 0 0 0 0 0 0 0 1027 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 thermsums
S 1027 11 0 0 0 8 986 582 7887 40800010 801000 0 192 0 0 1026 1026 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 analysis_m$12
S 1028 23 5 0 0 0 1030 582 7539 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 set_thermodynamic_sums
S 1029 1 3 0 0 261 1 1028 7901 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ths
S 1030 14 5 0 0 0 1 1028 7539 0 400000 0 0 0 48 1 0 0 0 0 0 0 0 0 0 0 0 0 82 0 582 0 0 0 0 set_thermodynamic_sums
F 1030 1 1029
S 1031 23 5 0 0 8 1032 582 7562 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_thermodynamic_sums
S 1032 14 5 0 0 261 1 1031 7562 4 400000 0 0 0 50 0 0 0 1033 0 0 0 0 0 0 0 0 0 86 0 582 0 0 0 0 get_thermodynamic_sums
F 1032 0
S 1033 1 3 0 0 261 1 1031 7562 14 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_thermodynamic_sums
S 1034 23 5 0 0 0 1035 582 7513 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 update_thermodynamic_sums
S 1035 14 5 0 0 0 1 1034 7513 0 400000 0 0 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 99 0 582 0 0 0 0 update_thermodynamic_sums
F 1035 0
S 1036 23 5 0 0 0 1037 582 7905 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bin
S 1037 14 5 0 0 0 1 1036 7905 10 400000 0 0 0 52 0 0 0 0 0 0 0 0 0 0 0 0 0 124 0 582 0 0 0 0 bin
F 1037 0
S 1038 23 5 0 0 0 1040 582 7585 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 runavs
S 1039 1 3 0 0 6 1 1038 7909 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n1
S 1040 14 5 0 0 0 1 1038 7585 0 400000 0 0 0 53 1 0 0 0 0 0 0 0 0 0 0 0 0 193 0 582 0 0 0 0 runavs
F 1040 1 1039
S 1041 23 5 0 0 0 1043 582 7596 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 auto
S 1042 1 3 0 0 6 1 1041 7912 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ifft
S 1043 14 5 0 0 0 1 1041 7596 0 400000 0 0 0 55 1 0 0 0 0 0 0 0 0 0 0 0 0 308 0 582 0 0 0 0 auto
F 1043 1 1042
S 1044 23 5 0 0 0 1047 582 7592 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rdf
S 1045 1 3 0 0 6 1 1044 7917 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 maxbin
S 1046 1 3 0 0 6 1 1044 7924 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 binsperangstrom
S 1047 14 5 0 0 0 1 1044 7592 0 400000 0 0 0 57 2 0 0 0 0 0 0 0 0 0 0 0 0 419 0 582 0 0 0 0 rdf
F 1047 2 1045 1046
A 15 2 0 0 0 6 592 0 0 0 15 0 0 0 0 0 0 0 0 0
A 25 2 0 0 0 6 595 0 0 0 25 0 0 0 0 0 0 0 0 0
A 52 2 0 0 0 6 637 0 0 0 52 0 0 0 0 0 0 0 0 0
A 56 2 0 0 0 6 632 0 0 0 56 0 0 0 0 0 0 0 0 0
A 57 2 0 0 0 6 633 0 0 0 57 0 0 0 0 0 0 0 0 0
A 58 2 0 0 0 9 577 0 0 0 58 0 0 0 0 0 0 0 0 0
A 59 2 0 0 0 8 573 0 0 0 59 0 0 0 0 0 0 0 0 0
A 60 2 0 0 0 6 634 0 0 0 60 0 0 0 0 0 0 0 0 0
A 61 2 0 0 0 16 635 0 0 0 61 0 0 0 0 0 0 0 0 0
A 62 2 0 0 0 6 636 0 0 0 62 0 0 0 0 0 0 0 0 0
Z
T 639 33 0 3 0 0
A 646 6 0 0 1 56 1
A 647 6 0 0 1 56 1
A 648 6 0 0 1 56 1
A 649 6 0 0 1 57 1
A 650 6 0 0 1 2 1
A 651 6 0 0 1 25 1
A 652 6 0 0 1 25 1
A 653 6 0 0 1 25 1
A 654 6 0 0 1 2 1
A 655 6 0 0 1 2 1
A 656 9 0 0 1 58 1
A 657 9 0 0 1 58 1
A 662 9 0 0 1 58 1
A 663 9 0 0 1 58 1
A 665 8 0 0 1 59 1
A 666 6 0 0 1 3 1
A 667 6 0 0 1 2 1
A 668 6 0 0 1 60 1
A 670 6 0 0 1 56 1
A 671 6 0 0 1 2 1
A 672 8 0 0 1 59 1
A 673 16 0 0 1 61 1
A 674 6 0 0 1 62 1
A 675 6 0 0 1 2 1
A 676 16 0 0 1 61 1
A 679 6 0 0 1 2 1
A 680 6 0 0 1 2 1
A 681 6 0 0 1 2 1
A 682 6 0 0 1 2 1
A 683 6 0 0 1 2 1
A 684 6 0 0 1 2 1
A 686 16 0 0 1 61 0
T 1009 261 0 3 0 0
A 1010 9 0 0 1 58 1
A 1011 9 0 0 1 58 1
A 1012 9 0 0 1 58 1
A 1013 9 0 0 1 58 1
A 1014 9 0 0 1 58 1
A 1015 9 0 0 1 58 1
A 1016 9 0 0 1 58 1
A 1017 9 0 0 1 58 1
A 1018 9 0 0 1 58 1
A 1019 9 0 0 1 58 1
A 1020 9 0 0 1 58 1
A 1021 9 0 0 1 58 1
A 1022 9 0 0 1 58 1
A 1023 9 0 0 1 58 1
A 1024 9 0 0 1 58 1
R 1025 267 0 0
A 0 9 0 52 1 58 0
Z
