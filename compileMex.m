function compileMex()
cd ..
mex private/edgesDetect_3D_Mex.cpp -outdir private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex private/houghVotes_3D_Mex.cpp  -outdir private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex private/edgesNmsMex.cpp        -outdir private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex toolbox/classify/private/regForestFindThr.cpp -outdir toolbox/classify/private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"