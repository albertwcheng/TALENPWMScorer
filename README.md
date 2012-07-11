TALENPWMScorer
==============

Dependencies
============
 * g++
 * python

Build
=====

    git clone https://github.com/andynu/TALENPWMScorer.git
    git submodule init
    git submodule update
    ./make.sh

Run
===

    Usage: ./findTALENPairsInGenome.py seqDir mode1 motifDef1 modeParam1 topN1 minSpacer maxSpacer mode2 motifDef2 modeParam2 topN2 outDir
    seqDir 
    mode1 
      "consensus"
        provided desired bind sequence
      "PWM"
    motifDef1 
      if mode1 "consensus"
        provided desired bind sequence
      if mode1 "PWM"
      filename that provides the position weight matrix
    modeParam1 
      if mode1 "consensus"
        number of mismatches 
      if mode1 "PWM"
        minimum score
    topN1 
      number of sites to output
    minSpacer 
      number
    maxSpacer 
      number
    mode2
      "consensus"
        provided desired bind sequence
      "PWM"
    motifDef2
      if mode2 "consensus"
        provided desired bind sequence
      if mode2 "PWM"
      filename that provides the position weight matrix
    modeParam2 
      if mode2 "consensus"
        number of mismatches 
      if mode2 "PWM"
        minimum score
    topN2g 
      number of sites to output
    minSpacer 
      number
    maxSpacer 
      number
    outDir
