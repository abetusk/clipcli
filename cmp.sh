#!/bin/bash

g++ -g -I. clipper.cpp sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc clipcli.cpp -o clipcli
