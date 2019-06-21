#!/usr/bin/env python

from brainiac import ProcessChecker
from brainiac import DicomConverter
from brainiac import NiftiRotator
from brainiac import ScanEprimeConverter
from brainiac import ScanEprimeDivider
from brainiac import OnsetCreator
from brainiac import BetBrainExtractor
from brainiac import BSEBrainExtractor
from brainiac import BiasCorrector
from brainiac import VolumeTrimmer
from brainiac import MotionEvaluator
from brainiac import JsonCreator


import sys
from datetime import datetime

startTime = datetime.now()

study_name = sys.argv[1]

subject_number = sys.argv[2]
subject_number = str(subject_number)

process_checker = ProcessChecker(study_name, subject_number)
process_checker.process()

nifti_converter = DicomConverter(study_name, subject_number)
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

json_creator = JsonCreator(study_name, subject_number)
json_creator.process()

print (datetime.now() - startTime)

