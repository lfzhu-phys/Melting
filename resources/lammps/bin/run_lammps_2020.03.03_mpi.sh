#!/bin/bash
mpiexec -n $1 --oversubscribe lmp_mpi -in control.inp;
