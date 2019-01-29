#!/usr/bin/env python

from braniac import ProcessChecker
from braniac import NiftiConverter
from braniac import NiftiRotator
from braniac import ScanEprimeConverter
from braniac import ScanEprimeDivider
from braniac import OnsetCreator
from braniac import BetBrainExtractor
from braniac import BSEBrainExtractor
from braniac import BiasCorrector
from braniac import VolumeTrimmer
from braniac import MotionEvaluator

import sys
from datetime import datetime

startTime = datetime.now()

study_name = sys.argv[1]

subject_number = sys.argv[2]
subject_number = str(subject_number)

process_checker = ProcessChecker(study_name, subject_number)
process_checker.process()

nifti_converter = NiftiConverter(study_name, subject_number)
nifti_converter.process()

rotator = NiftiRotator(study_name, subject_number)
rotator.process()

eprime_converter = ScanEprimeConverter(study_name, subject_number)
eprime_converter.process()

eprime_divider = ScanEprimeDivider(study_name, subject_number)
eprime_divider.process()

onset_creator = OnsetCreator(study_name, subject_number)
onset_creator.process()

bet_brain_extractor = BetBrainExtractor(study_name, subject_number)
bet_brain_extractor.process()

bse_brain_extractor = BSEBrainExtractor(study_name, subject_number)
bse_brain_extractor.process()

bias_corrector = BiasCorrector(study_name, subject_number)
bias_corrector.process()

volume_trimmer = VolumeTrimmer(study_name, subject_number)
volume_trimmer.process()

motion_evaluator = MotionEvaluator(study_name, subject_number)
motion_evaluator.process()

print (datetime.now() - startTime)

