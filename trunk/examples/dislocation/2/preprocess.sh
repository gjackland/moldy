gfortran -O4 mod_disloc.f90 disloc.f90 -o disloc.exe
./disloc.exe system.in system.in.new ../0/system.in
rm system.in
mv system.in.new system.in