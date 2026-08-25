#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : IndexMAF
# @created     : Saturday Feb 08, 2025 15:28:26 CST
#
# @description : 
######################################################################

import sys
from Bio.AlignIO import MafIO
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    # sys.argv = ["",  "multiple_sequence_alignment_analysis/maf/chr14_KI270723v1_random.maf","multiple_sequence_alignment_analysis/maf/chr14_KI270723v1_random.mafindex", "hg38.chr14_KI270723v1_random"]
    sys.argv = ["",  "multiple_sequence_alignment_analysis/maf/chr14_KI270723v1_random.maf","multiple_sequence_alignment_analysis/maf/chr14_KI270723v1_random.mafindex", "hg38.chr14_KI270723v1_random"]

_, MAF, MAFindex, seqname = sys.argv

idx = MafIO.MafIndex(MAFindex ,MAF, seqname)
