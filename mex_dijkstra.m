% Create matlab executable of dijkstra.cpp

% bsub -N matlab -nodisplay -nojvm -singleCompThread -r mex_dijkstra


mex -largeArrayDims dijkstra.cpp