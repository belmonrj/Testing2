#!/bin/sh

g++ -o test_prod.exe test_prod.C ntp_event.C -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

