V24 timers_m
12 timers_m.F90 S582 0
08/28/2009  14:35:56
use params_m private
use params_m private
enduse
D 33 24 600 440 595 7
S 582 24 0 0 0 8 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 26 0 0 0 0 0 0 timers_m
S 586 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 588 3 0 0 0 6 1 1 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 589 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 590 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 591 3 0 0 0 16 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 592 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 595 25 1 params_m simparameters
R 600 5 6 params_m title1 simparameters
R 601 5 7 params_m title2 simparameters
R 602 5 8 params_m ivol simparameters
R 603 5 9 params_m iquen simparameters
R 604 5 10 params_m iprint simparameters
R 605 5 11 params_m nnbrs simparameters
R 606 5 12 params_m iverlet simparameters
R 607 5 13 params_m nlcx simparameters
R 608 5 14 params_m nlcy simparameters
R 609 5 15 params_m nlcz simparameters
R 610 5 16 params_m nm simparameters
R 611 5 17 params_m nspec simparameters
R 612 5 18 params_m dsp simparameters
R 613 5 19 params_m rpad simparameters
R 614 5 20 params_m rcut simparameters
R 615 5 21 params_m rnear simparameters
R 616 5 22 params_m deltat simparameters
R 617 5 23 params_m temprq simparameters
R 618 5 24 params_m rqke simparameters
R 619 5 25 params_m press simparameters
R 620 5 26 params_m nloops simparameters
R 621 5 27 params_m nsteps simparameters
R 622 5 28 params_m nprint simparameters
R 623 5 29 params_m ntcm simparameters
R 624 5 30 params_m nchkpt simparameters
R 625 5 31 params_m restart simparameters
R 626 5 32 params_m nose simparameters
R 627 5 33 params_m alternate_quench_md simparameters
R 628 5 34 params_m nout simparameters
R 629 5 35 params_m ntape simparameters
R 630 5 36 params_m dumpx1 simparameters
R 631 5 37 params_m tjob simparameters
R 632 5 38 params_m tfinalise simparameters
R 633 5 39 params_m prevsteps simparameters
R 634 5 40 params_m currentstep simparameters
R 635 5 41 params_m laststep simparameters
R 636 5 42 params_m lastprint simparameters
R 637 5 43 params_m lastchkpt simparameters
R 638 5 44 params_m ntc simparameters
R 639 5 45 params_m strx simparameters
R 640 5 46 params_m uselookup simparameters
R 641 5 47 params_m boxtem simparameters
R 642 5 48 params_m bdel2 simparameters
R 643 5 49 params_m bmass simparameters
S 707 27 0 0 0 8 710 582 5454 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 clock
S 708 6 4 0 0 9 1 582 5460 80000c 0 0 0 0 0 0 0 0 0 709 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 omptime
S 709 11 0 0 0 8 653 582 5471 40800000 801000 0 8 0 0 708 708 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 timers_m$10
S 710 23 5 0 0 8 711 582 5454 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 clock
S 711 14 5 0 0 9 1 710 5454 4 400000 0 0 0 11 0 0 0 712 0 0 0 0 0 0 0 0 0 42 0 582 0 0 0 0 clock
F 711 0
S 712 1 3 0 0 9 1 710 5454 14 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 clock
A 17 2 0 0 0 6 586 0 0 0 17 0 0 0 0 0 0 0 0 0
A 18 2 0 0 0 6 588 0 0 0 18 0 0 0 0 0 0 0 0 0
A 19 2 0 0 0 6 589 0 0 0 19 0 0 0 0 0 0 0 0 0
A 20 2 0 0 0 9 577 0 0 0 20 0 0 0 0 0 0 0 0 0
A 21 2 0 0 0 8 573 0 0 0 21 0 0 0 0 0 0 0 0 0
A 22 2 0 0 0 6 590 0 0 0 22 0 0 0 0 0 0 0 0 0
A 23 2 0 0 0 16 591 0 0 0 23 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 6 592 0 0 0 24 0 0 0 0 0 0 0 0 0
Z
T 595 33 0 3 0 0
A 602 6 0 0 1 18 1
A 603 6 0 0 1 18 1
A 604 6 0 0 1 18 1
A 605 6 0 0 1 19 1
A 606 6 0 0 1 2 1
A 607 6 0 0 1 17 1
A 608 6 0 0 1 17 1
A 609 6 0 0 1 17 1
A 610 6 0 0 1 2 1
A 611 6 0 0 1 2 1
A 612 9 0 0 1 20 1
A 613 9 0 0 1 20 1
A 619 8 0 0 1 21 1
A 620 6 0 0 1 3 1
A 621 6 0 0 1 2 1
A 622 6 0 0 1 22 1
A 624 6 0 0 1 18 1
A 625 6 0 0 1 2 1
A 626 8 0 0 1 21 1
A 627 16 0 0 1 23 1
A 628 6 0 0 1 24 1
A 629 6 0 0 1 2 1
A 630 16 0 0 1 23 1
A 633 6 0 0 1 2 1
A 634 6 0 0 1 2 1
A 635 6 0 0 1 2 1
A 636 6 0 0 1 2 1
A 637 6 0 0 1 2 1
A 638 6 0 0 1 2 1
A 640 16 0 0 1 23 0
Z
