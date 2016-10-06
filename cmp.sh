#!/bin/bash

#g++ -g -I. clipper.cpp sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc clipcli.cpp -o clipcli
g++ -std=c++11 -g -I. -I./lib clipper.cpp lib/delaunay.cpp lib/triangle.cpp clipcli.cpp -o clipcli
