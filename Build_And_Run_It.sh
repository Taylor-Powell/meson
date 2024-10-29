cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

mpiexec -n 1 build/meson data/b1.Ecm.dat