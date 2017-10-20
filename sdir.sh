#!/usr/bin/env bash

function sdir() {
  if [ -z "$1" ]
  then
    echo "Please provide slurm jobid"
  else
    DIR=`scontrol show jobid -dd $1 | grep WorkDir | sed 's/   WorkDir=//g'`
    cd $DIR
    ls
  fi
}
