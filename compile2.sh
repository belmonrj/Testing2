#!/bin/sh

g++ -o testing2.exe testing2.C ntp_event.C -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

