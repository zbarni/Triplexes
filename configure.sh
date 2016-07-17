#!/bin/bash
# This file needs to be sourced instead of executed.

echo "Setting TRIPLEXATOR_HOME env. variable"
export TRIPLEXATOR_HOME=$PWD
echo "Setting TRIPLEXATOR_LIBRARY env. variable"
export TRIPLEXATOR_LIBRARY=$TRIPLEXATOR_HOME/triplexator/lib/libtriplexator.so

