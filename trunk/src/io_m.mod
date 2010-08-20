V24 io_m
8 io_m.F90 S582 0
04/26/2010  18:05:26
use params_m private
use params_m private
enduse
D 36 24 658 456 653 7
S 582 24 0 0 0 6 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 26 0 0 0 0 0 0 io_m
S 594 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 646 3 0 0 0 6 1 1 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 647 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 648 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 649 3 0 0 0 16 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 650 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 653 25 1 params_m simparameters
R 658 5 6 params_m title1 simparameters
R 659 5 7 params_m title2 simparameters
R 660 5 8 params_m ivol simparameters
R 661 5 9 params_m iquen simparameters
R 662 5 10 params_m iprint simparameters
R 663 5 11 params_m nnbrs simparameters
R 664 5 12 params_m iverlet simparameters
R 665 5 13 params_m nlcx simparameters
R 666 5 14 params_m nlcy simparameters
R 667 5 15 params_m nlcz simparameters
R 668 5 16 params_m nm simparameters
R 669 5 17 params_m nspec simparameters
R 670 5 18 params_m dsp simparameters
R 671 5 19 params_m rpad simparameters
R 672 5 20 params_m rcut simparameters
R 673 5 21 params_m rnear simparameters
R 674 5 22 params_m deltat simparameters
R 675 5 23 params_m temprq simparameters
R 676 5 24 params_m tempsp simparameters
R 677 5 25 params_m zlayer simparameters
R 678 5 26 params_m rqke simparameters
R 679 5 27 params_m press simparameters
R 680 5 28 params_m nloops simparameters
R 681 5 29 params_m nsteps simparameters
R 682 5 30 params_m nprint simparameters
R 683 5 31 params_m ntcm simparameters
R 684 5 32 params_m nchkpt simparameters
R 685 5 33 params_m restart simparameters
R 686 5 34 params_m nose simparameters
R 687 5 35 params_m alternate_quench_md simparameters
R 688 5 36 params_m nout simparameters
R 689 5 37 params_m ntape simparameters
R 690 5 38 params_m dumpx1 simparameters
R 691 5 39 params_m tjob simparameters
R 692 5 40 params_m tfinalise simparameters
R 693 5 41 params_m prevsteps simparameters
R 694 5 42 params_m currentstep simparameters
R 695 5 43 params_m laststep simparameters
R 696 5 44 params_m lastprint simparameters
R 697 5 45 params_m lastchkpt simparameters
R 698 5 46 params_m ntc simparameters
R 699 5 47 params_m strx simparameters
R 700 5 48 params_m uselookup simparameters
R 701 5 49 params_m boxtem simparameters
R 702 5 50 params_m bdel2 simparameters
R 703 5 51 params_m bmass simparameters
S 956 27 0 0 0 8 964 582 7155 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 read_system
S 957 27 0 0 0 8 966 582 7167 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 write_textout_header
S 958 27 0 0 0 8 969 582 7188 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 write_textout
S 959 6 4 0 0 16 960 582 7202 14 0 0 0 0 0 0 0 0 0 962 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 exists
S 960 6 4 0 0 16 1 582 7209 14 0 0 4 0 0 0 0 0 0 962 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 opened
S 961 6 4 0 0 36 1 582 5500 80003c 0 0 0 0 0 0 0 0 0 963 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 simparam
S 962 11 0 0 0 6 931 582 7216 40800010 801000 0 8 0 0 959 960 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 io_m$4
S 963 11 0 0 0 6 962 582 7223 40800010 801000 0 456 0 0 961 961 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 io_m$12
S 964 23 5 0 0 0 965 582 7155 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 read_system
S 965 14 5 0 0 0 1 964 7155 0 400000 0 0 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 60 0 582 0 0 0 0 read_system
F 965 0
S 966 23 5 0 0 0 968 582 7167 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 write_textout_header
S 967 1 3 1 0 9 1 966 7231 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 local_vol
S 968 14 5 0 0 0 1 966 7167 0 400000 0 0 0 46 1 0 0 0 0 0 0 0 0 0 0 0 0 210 0 582 0 0 0 0 write_textout_header
F 968 1 967
S 969 23 5 0 0 0 970 582 7188 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 write_textout
S 970 14 5 0 0 0 1 969 7188 0 400000 0 0 0 48 0 0 0 0 0 0 0 0 0 0 0 0 0 302 0 582 0 0 0 0 write_textout
F 970 0
A 25 2 0 0 0 6 594 0 0 0 25 0 0 0 0 0 0 0 0 0
A 60 2 0 0 0 6 646 0 0 0 60 0 0 0 0 0 0 0 0 0
A 61 2 0 0 0 6 647 0 0 0 61 0 0 0 0 0 0 0 0 0
A 62 2 0 0 0 9 577 0 0 0 62 0 0 0 0 0 0 0 0 0
A 63 2 0 0 0 8 573 0 0 0 63 0 0 0 0 0 0 0 0 0
A 64 2 0 0 0 6 648 0 0 0 64 0 0 0 0 0 0 0 0 0
A 65 2 0 0 0 16 649 0 0 0 65 0 0 0 0 0 0 0 0 0
A 66 2 0 0 0 6 650 0 0 0 66 0 0 0 0 0 0 0 0 0
Z
T 653 36 0 3 0 0
A 660 6 0 0 1 60 1
A 661 6 0 0 1 60 1
A 662 6 0 0 1 60 1
A 663 6 0 0 1 61 1
A 664 6 0 0 1 2 1
A 665 6 0 0 1 25 1
A 666 6 0 0 1 25 1
A 667 6 0 0 1 25 1
A 668 6 0 0 1 2 1
A 669 6 0 0 1 2 1
A 670 9 0 0 1 62 1
A 671 9 0 0 1 62 1
A 676 9 0 0 1 62 1
A 677 9 0 0 1 62 1
A 679 8 0 0 1 63 1
A 680 6 0 0 1 3 1
A 681 6 0 0 1 2 1
A 682 6 0 0 1 64 1
A 684 6 0 0 1 60 1
A 685 6 0 0 1 2 1
A 686 8 0 0 1 63 1
A 687 16 0 0 1 65 1
A 688 6 0 0 1 66 1
A 689 6 0 0 1 2 1
A 690 16 0 0 1 65 1
A 693 6 0 0 1 2 1
A 694 6 0 0 1 2 1
A 695 6 0 0 1 2 1
A 696 6 0 0 1 2 1
A 697 6 0 0 1 2 1
A 698 6 0 0 1 2 1
A 700 16 0 0 1 65 0
Z
