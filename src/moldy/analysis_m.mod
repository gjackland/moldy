V24 analysis_m
14 analysis_m.F90 S582 0
09/22/2009  09:28:48
use params_m private
use zirconium private
use params_m private
use zirconium private
enduse
D 33 24 645 448 640 7
D 260 21 6 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 263 21 6 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 326 24 1161 192 1160 7
D 332 21 9 2 15 52 0 0 0 0 0
 0 25 3 3 25 25
 0 25 25 3 25 25
S 582 24 0 0 0 8 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 27 0 0 0 0 0 0 analysis_m
S 593 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 596 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 633 3 0 0 0 6 1 1 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 634 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 635 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 636 3 0 0 0 16 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 637 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 638 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 640 25 1 params_m simparameters
R 645 5 6 params_m title1 simparameters
R 646 5 7 params_m title2 simparameters
R 647 5 8 params_m ivol simparameters
R 648 5 9 params_m iquen simparameters
R 649 5 10 params_m iprint simparameters
R 650 5 11 params_m nnbrs simparameters
R 651 5 12 params_m iverlet simparameters
R 652 5 13 params_m nlcx simparameters
R 653 5 14 params_m nlcy simparameters
R 654 5 15 params_m nlcz simparameters
R 655 5 16 params_m nm simparameters
R 656 5 17 params_m nspec simparameters
R 657 5 18 params_m dsp simparameters
R 658 5 19 params_m rpad simparameters
R 659 5 20 params_m rcut simparameters
R 660 5 21 params_m rnear simparameters
R 661 5 22 params_m deltat simparameters
R 662 5 23 params_m temprq simparameters
R 663 5 24 params_m tempsp simparameters
R 664 5 25 params_m rqke simparameters
R 665 5 26 params_m press simparameters
R 666 5 27 params_m nloops simparameters
R 667 5 28 params_m nsteps simparameters
R 668 5 29 params_m nprint simparameters
R 669 5 30 params_m ntcm simparameters
R 670 5 31 params_m nchkpt simparameters
R 671 5 32 params_m restart simparameters
R 672 5 33 params_m nose simparameters
R 673 5 34 params_m alternate_quench_md simparameters
R 674 5 35 params_m nout simparameters
R 675 5 36 params_m ntape simparameters
R 676 5 37 params_m dumpx1 simparameters
R 677 5 38 params_m tjob simparameters
R 678 5 39 params_m tfinalise simparameters
R 679 5 40 params_m prevsteps simparameters
R 680 5 41 params_m currentstep simparameters
R 681 5 42 params_m laststep simparameters
R 682 5 43 params_m lastprint simparameters
R 683 5 44 params_m lastchkpt simparameters
R 684 5 45 params_m ntc simparameters
R 685 5 46 params_m strx simparameters
R 686 5 47 params_m uselookup simparameters
R 687 5 48 params_m boxtem simparameters
R 688 5 49 params_m bdel2 simparameters
R 689 5 50 params_m bmass simparameters
S 1006 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 1022 7 16 zirconium supported_atomic_numbers$ac
S 1154 27 0 0 0 8 1185 582 8392 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 update_thermodynamic_sums
S 1155 27 0 0 0 8 1179 582 8418 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 set_thermodynamic_sums
S 1156 27 0 0 0 8 1182 582 8441 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 get_thermodynamic_sums
S 1157 27 0 0 0 8 1189 582 8464 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 runavs
S 1158 27 0 0 0 8 1195 582 8471 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 rdf
S 1159 27 0 0 0 8 1192 582 8475 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 auto
S 1160 25 0 0 0 326 1 582 8480 c 800000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 thermodynamic_sums
S 1161 5 0 0 0 9 1162 582 8499 800004 0 0 0 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1 1161 0 582 0 0 0 0 spe
S 1162 5 0 0 0 9 1163 582 8514 800004 0 0 8 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1161 1162 0 582 0 0 0 0 ske
S 1163 5 0 0 0 9 1164 582 8529 800004 0 0 16 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1162 1163 0 582 0 0 0 0 ste
S 1164 5 0 0 0 9 1165 582 8544 800004 0 0 24 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1163 1164 0 582 0 0 0 0 sh
S 1165 5 0 0 0 9 1166 582 8558 800004 0 0 32 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1164 1165 0 582 0 0 0 0 sth
S 1166 5 0 0 0 9 1167 582 8573 800004 0 0 40 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1165 1166 0 582 0 0 0 0 sf2
S 1167 5 0 0 0 9 1168 582 8588 800004 0 0 48 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1166 1167 0 582 0 0 0 0 svol
S 1168 5 0 0 0 9 1169 582 8604 800004 0 0 56 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1167 1168 0 582 0 0 0 0 spesq
S 1169 5 0 0 0 9 1170 582 8621 800004 0 0 64 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1168 1169 0 582 0 0 0 0 skesq
S 1170 5 0 0 0 9 1171 582 8638 800004 0 0 72 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1169 1170 0 582 0 0 0 0 stesq
S 1171 5 0 0 0 9 1172 582 8655 800004 0 0 80 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1170 1171 0 582 0 0 0 0 shsq
S 1172 5 0 0 0 9 1173 582 8671 800004 0 0 88 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1171 1172 0 582 0 0 0 0 sthsq
S 1173 5 0 0 0 9 1174 582 8688 800004 0 0 96 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1172 1173 0 582 0 0 0 0 sf2sq
S 1174 5 0 0 0 9 1175 582 8705 800004 0 0 104 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1173 1174 0 582 0 0 0 0 svolsq
S 1175 5 0 0 0 9 1176 582 8723 800004 0 0 112 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1174 1175 0 582 0 0 0 0 svolpe
S 1176 5 0 0 0 332 1 582 8741 800004 0 0 120 0 0 326 0 0 0 0 0 0 0 0 0 0 0 1175 1176 0 582 0 0 0 0 sb0
S 1177 6 4 0 0 326 1 582 8756 80003c 0 0 0 0 0 0 0 0 0 1178 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 thermsums
S 1178 11 0 0 0 8 1131 582 8766 40800010 801000 0 192 0 0 1177 1177 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 analysis_m$12
S 1179 23 5 0 0 0 1181 582 8418 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 set_thermodynamic_sums
S 1180 1 3 0 0 326 1 1179 8780 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ths
S 1181 14 5 0 0 0 1 1179 8418 0 400000 0 0 0 92 1 0 0 0 0 0 0 0 0 0 0 0 0 82 0 582 0 0 0 0 set_thermodynamic_sums
F 1181 1 1180
S 1182 23 5 0 0 8 1183 582 8441 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_thermodynamic_sums
S 1183 14 5 0 0 326 1 1182 8441 4 400000 0 0 0 94 0 0 0 1184 0 0 0 0 0 0 0 0 0 86 0 582 0 0 0 0 get_thermodynamic_sums
F 1183 0
S 1184 1 3 0 0 326 1 1182 8441 14 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_thermodynamic_sums
S 1185 23 5 0 0 0 1186 582 8392 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 update_thermodynamic_sums
S 1186 14 5 0 0 0 1 1185 8392 0 400000 0 0 0 95 0 0 0 0 0 0 0 0 0 0 0 0 0 99 0 582 0 0 0 0 update_thermodynamic_sums
F 1186 0
S 1187 23 5 0 0 0 1188 582 8784 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bin
S 1188 14 5 0 0 0 1 1187 8784 10 400000 0 0 0 96 0 0 0 0 0 0 0 0 0 0 0 0 0 124 0 582 0 0 0 0 bin
F 1188 0
S 1189 23 5 0 0 0 1191 582 8464 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 runavs
S 1190 1 3 0 0 6 1 1189 8788 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n1
S 1191 14 5 0 0 0 1 1189 8464 0 400000 0 0 0 97 1 0 0 0 0 0 0 0 0 0 0 0 0 193 0 582 0 0 0 0 runavs
F 1191 1 1190
S 1192 23 5 0 0 0 1194 582 8475 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 auto
S 1193 1 3 0 0 6 1 1192 8791 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ifft
S 1194 14 5 0 0 0 1 1192 8475 0 400000 0 0 0 99 1 0 0 0 0 0 0 0 0 0 0 0 0 293 0 582 0 0 0 0 auto
F 1194 1 1193
S 1195 23 5 0 0 0 1198 582 8471 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rdf
S 1196 1 3 0 0 6 1 1195 8796 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 maxbin
S 1197 1 3 0 0 6 1 1195 8803 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 binsperangstrom
S 1198 14 5 0 0 0 1 1195 8471 0 400000 0 0 0 101 2 0 0 0 0 0 0 0 0 0 0 0 0 404 0 582 0 0 0 0 rdf
F 1198 2 1196 1197
A 15 2 0 0 0 6 593 0 0 0 15 0 0 0 0 0 0 0 0 0
A 25 2 0 0 0 6 596 0 0 0 25 0 0 0 0 0 0 0 0 0
A 52 2 0 0 0 6 638 0 0 0 52 0 0 0 0 0 0 0 0 0
A 56 2 0 0 0 6 633 0 0 0 56 0 0 0 0 0 0 0 0 0
A 57 2 0 0 0 6 634 0 0 0 57 0 0 0 0 0 0 0 0 0
A 58 2 0 0 0 9 577 0 0 0 58 0 0 0 0 0 0 0 0 0
A 59 2 0 0 0 8 573 0 0 0 59 0 0 0 0 0 0 0 0 0
A 60 2 0 0 0 6 635 0 0 0 60 0 0 0 0 0 0 0 0 0
A 61 2 0 0 0 16 636 0 0 0 61 0 0 0 0 0 0 0 0 0
A 62 2 0 0 0 6 637 0 0 0 62 0 0 0 0 0 0 0 0 0
A 282 2 0 0 21 6 1006 0 0 0 282 0 0 0 0 0 0 0 0 0
A 290 1 0 5 0 260 1022 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 67 1 1
V 290 260 7 0
R 0 263 0 0
A 0 6 0 0 1 282 0
T 640 33 0 3 0 0
A 647 6 0 0 1 56 1
A 648 6 0 0 1 56 1
A 649 6 0 0 1 56 1
A 650 6 0 0 1 57 1
A 651 6 0 0 1 2 1
A 652 6 0 0 1 25 1
A 653 6 0 0 1 25 1
A 654 6 0 0 1 25 1
A 655 6 0 0 1 2 1
A 656 6 0 0 1 2 1
A 657 9 0 0 1 58 1
A 658 9 0 0 1 58 1
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
T 1160 326 0 3 0 0
A 1161 9 0 0 1 58 1
A 1162 9 0 0 1 58 1
A 1163 9 0 0 1 58 1
A 1164 9 0 0 1 58 1
A 1165 9 0 0 1 58 1
A 1166 9 0 0 1 58 1
A 1167 9 0 0 1 58 1
A 1168 9 0 0 1 58 1
A 1169 9 0 0 1 58 1
A 1170 9 0 0 1 58 1
A 1171 9 0 0 1 58 1
A 1172 9 0 0 1 58 1
A 1173 9 0 0 1 58 1
A 1174 9 0 0 1 58 1
A 1175 9 0 0 1 58 1
R 1176 332 0 0
A 0 9 0 52 1 58 0
Z
