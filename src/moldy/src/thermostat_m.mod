V24 thermostat_m
16 thermostat_m.F90 S582 0
08/28/2009  14:35:55
use params_m public 0 direct
use system_m public 0 direct
use constants_m public 0 direct
use random_m public 0 direct
use particles_m public 0 direct
use parinellorahman_m public 0 direct
enduse
D 33 24 604 440 599 7
S 582 24 0 0 0 8 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 34 0 0 0 0 0 0 thermostat_m
S 590 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 592 3 0 0 0 6 1 1 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 593 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 594 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 595 3 0 0 0 16 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 596 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 599 25 1 params_m simparameters
R 604 5 6 params_m title1 simparameters
R 605 5 7 params_m title2 simparameters
R 606 5 8 params_m ivol simparameters
R 607 5 9 params_m iquen simparameters
R 608 5 10 params_m iprint simparameters
R 609 5 11 params_m nnbrs simparameters
R 610 5 12 params_m iverlet simparameters
R 611 5 13 params_m nlcx simparameters
R 612 5 14 params_m nlcy simparameters
R 613 5 15 params_m nlcz simparameters
R 614 5 16 params_m nm simparameters
R 615 5 17 params_m nspec simparameters
R 616 5 18 params_m dsp simparameters
R 617 5 19 params_m rpad simparameters
R 618 5 20 params_m rcut simparameters
R 619 5 21 params_m rnear simparameters
R 620 5 22 params_m deltat simparameters
R 621 5 23 params_m temprq simparameters
R 622 5 24 params_m rqke simparameters
R 623 5 25 params_m press simparameters
R 624 5 26 params_m nloops simparameters
R 625 5 27 params_m nsteps simparameters
R 626 5 28 params_m nprint simparameters
R 627 5 29 params_m ntcm simparameters
R 628 5 30 params_m nchkpt simparameters
R 629 5 31 params_m restart simparameters
R 630 5 32 params_m nose simparameters
R 631 5 33 params_m alternate_quench_md simparameters
R 632 5 34 params_m nout simparameters
R 633 5 35 params_m ntape simparameters
R 634 5 36 params_m dumpx1 simparameters
R 635 5 37 params_m tjob simparameters
R 636 5 38 params_m tfinalise simparameters
R 637 5 39 params_m prevsteps simparameters
R 638 5 40 params_m currentstep simparameters
R 639 5 41 params_m laststep simparameters
R 640 5 42 params_m lastprint simparameters
R 641 5 43 params_m lastchkpt simparameters
R 642 5 44 params_m ntc simparameters
R 643 5 45 params_m strx simparameters
R 644 5 46 params_m uselookup simparameters
R 645 5 47 params_m boxtem simparameters
R 646 5 48 params_m bdel2 simparameters
R 647 5 49 params_m bmass simparameters
S 923 27 0 0 0 6 929 582 6984 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 init_thermostat_m
S 924 27 0 0 0 8 931 582 7002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 cleanup_thermostat_m
S 925 6 4 0 0 9 1 582 7023 80001c 0 0 0 0 0 0 0 0 0 927 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 snhv
S 926 6 4 0 0 9 1 582 7034 14 0 0 0 0 0 0 0 0 0 928 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 temp
S 927 11 0 0 0 8 903 582 7039 40800010 801000 0 8 0 0 925 925 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 thermostat_m$14
S 928 11 0 0 0 8 927 582 7055 40800010 801000 0 8 0 0 926 926 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 thermostat_m$6
S 929 23 5 0 0 0 930 582 6984 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 init_thermostat_m
S 930 14 5 0 0 0 1 929 6984 0 400000 0 0 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 61 0 582 0 0 0 0 init_thermostat_m
F 930 0
S 931 23 5 0 0 0 932 582 7002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cleanup_thermostat_m
S 932 14 5 0 0 0 1 931 7002 0 400000 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 67 0 582 0 0 0 0 cleanup_thermostat_m
F 932 0
S 933 23 5 0 0 0 935 582 7070 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 set_temp
S 934 1 3 1 0 9 1 933 5257 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dum
S 935 14 5 0 0 0 1 933 7070 0 400000 0 0 0 33 1 0 0 0 0 0 0 0 0 0 0 0 0 76 0 582 0 0 0 0 set_temp
F 935 1 934
S 936 23 5 0 0 8 937 582 7079 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_temp
S 937 14 5 0 0 9 1 936 7079 1004 400000 0 0 0 35 0 0 0 938 0 0 0 0 0 0 0 0 0 80 0 582 0 0 0 0 get_temp
F 937 0
S 938 1 3 0 0 9 1 936 7079 4 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_temp
S 939 23 5 0 0 0 940 582 7088 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 kinten
S 940 14 5 0 0 0 1 939 7088 0 400000 0 0 0 36 0 0 0 0 0 0 0 0 0 0 0 0 0 93 0 582 0 0 0 0 kinten
F 940 0
S 941 23 5 0 0 0 942 582 7095 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setvel
S 942 14 5 0 0 0 1 941 7095 0 400000 0 0 0 37 0 0 0 0 0 0 0 0 0 0 0 0 0 159 0 582 0 0 0 0 setvel
F 942 0
S 943 23 5 0 0 0 944 582 7102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 scalev
S 944 14 5 0 0 0 1 943 7102 0 400000 0 0 0 38 0 0 0 0 0 0 0 0 0 0 0 0 0 284 0 582 0 0 0 0 scalev
F 944 0
A 17 2 0 0 0 6 590 0 0 0 17 0 0 0 0 0 0 0 0 0
A 18 2 0 0 0 6 592 0 0 0 18 0 0 0 0 0 0 0 0 0
A 19 2 0 0 0 6 593 0 0 0 19 0 0 0 0 0 0 0 0 0
A 20 2 0 0 0 9 577 0 0 0 20 0 0 0 0 0 0 0 0 0
A 21 2 0 0 0 8 573 0 0 0 21 0 0 0 0 0 0 0 0 0
A 22 2 0 0 0 6 594 0 0 0 22 0 0 0 0 0 0 0 0 0
A 23 2 0 0 0 16 595 0 0 0 23 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 6 596 0 0 0 24 0 0 0 0 0 0 0 0 0
Z
T 599 33 0 3 0 0
A 606 6 0 0 1 18 1
A 607 6 0 0 1 18 1
A 608 6 0 0 1 18 1
A 609 6 0 0 1 19 1
A 610 6 0 0 1 2 1
A 611 6 0 0 1 17 1
A 612 6 0 0 1 17 1
A 613 6 0 0 1 17 1
A 614 6 0 0 1 2 1
A 615 6 0 0 1 2 1
A 616 9 0 0 1 20 1
A 617 9 0 0 1 20 1
A 623 8 0 0 1 21 1
A 624 6 0 0 1 3 1
A 625 6 0 0 1 2 1
A 626 6 0 0 1 22 1
A 628 6 0 0 1 18 1
A 629 6 0 0 1 2 1
A 630 8 0 0 1 21 1
A 631 16 0 0 1 23 1
A 632 6 0 0 1 24 1
A 633 6 0 0 1 2 1
A 634 16 0 0 1 23 1
A 637 6 0 0 1 2 1
A 638 6 0 0 1 2 1
A 639 6 0 0 1 2 1
A 640 6 0 0 1 2 1
A 641 6 0 0 1 2 1
A 642 6 0 0 1 2 1
A 644 16 0 0 1 23 0
Z
