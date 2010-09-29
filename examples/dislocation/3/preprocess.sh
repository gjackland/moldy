gfortran -O4 freeze_surface.f90 -o freeze_surface.exe
./freeze_surface.exe system.in system.in.new
rm system.in
mv system.in.new system.in