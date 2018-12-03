#!/usr/bin/env python

"""testr.py: Updating the original nueroimaging processing pipeline"""

__author__ = "David Lee"
__credits__ = ["David Lee", "Micheal Kelly", "Jeanette Mumford", "Nate Vack"]
__version__ = "1.0"
__maintainer__ = "David Lee"
__email__ = "david.s.lee@wisc.edu"
__status__ = "Production"

import os
import glob
import csv
import pandas as pd
import shutil

# TO DO
# 1. add print statements in general
# 2. qa all the outputs. probably just compare and contrast with the old one

class MotionEvaluator:
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	def identify_qa_directory(self):
		qa_dir = f'{self.output_dir}{self.subject_dir}/func/QA'
		return qa_dir

	def create_qa_directory(self):
		qa_dir = self.identify_qa_directory()
		os.makedirs(qa_dir, exist_ok=True)

	def identify_html_file(self):
		qa_dir = self.identify_qa_directory()
		html_file = f'{qa_dir}/motion.html'
		return html_file

	def evaluate_task_fmri_motion(self):
		qa_dir = self.identify_qa_directory()
		html_file = self.identify_html_file()
		volume_trimmer = VolumeTrimmer(self.study_name, self.subject_number)
		for i in range(1,4):
			infile = volume_trimmer.identify_task_fmri_trimmed_file("EmotionRegulation", i)
			confound_outfile = f'{qa_dir}/0{i}_confound.txt'
			plot_outfile = f'{qa_dir}/fd_plot_0{i}.png'
			outlier_outfile = f'{qa_dir}/outlier_output_0{i}.txt'
			os.system(f'fsl_motion_outliers -i {infile} -o {confound_outfile} --fd --thresh=0.5 -p {plot_outfile} -v > {outlier_outfile}')
			if not os.path.isfile(confound_outfile):
				os.mknod(confound_outfile)
			os.system(f'cat {outlier_outfile} >> {html_file}')
			os.system(f"echo '<br><br> FD plot {self.subject_number} 0{i} <br><IMG BORDER=0 SRC={plot_outfile} WIDTH=100%></BODY></HTML> <p>==========================================================<p>' >> {html_file}")

	def evaluate_resting_fmri_motion(self):
		qa_dir = self.identify_qa_directory()
		html_file = self.identify_html_file()
		volume_trimmer = VolumeTrimmer(self.study_name, self.subject_number)
		infile = volume_trimmer.identify_resting_fmri_trimmed_file()
		confound_outfile = f'{qa_dir}/04_confound.txt'
		plot_outfile = f'{qa_dir}/fd_plot_04.png'
		outlier_outfile = f'{qa_dir}/outlier_output_04.txt'
		os.system(f'fsl_motion_outliers -i {infile} -o {confound_outfile} --fd --thresh=0.2 -p {plot_outfile} -v > {outlier_outfile}')
		if not os.path.isfile(confound_outfile):
				os.mknod(confound_outfile)
		os.system(f'cat {outlier_outfile} >> {html_file}')
		os.system(f"echo '<br><br> FD plot {self.subject_number} 04 <br><IMG BORDER=0 SRC={plot_outfile} WIDTH=100%></BODY></HTML> <p>==========================================================<p>' >> {html_file}")


	def process(self):
		self.create_qa_directory()
		self.evaluate_task_fmri_motion()
		self.evaluate_resting_fmri_motion()

class VolumeTrimmer:

	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	def identify_task_fmri_input_file(self, task_type, run_number):
		infile = f'{self.input_dir}{self.subject_dir}/func/{self.subject_dir}_task-{task_type}_run-0{run_number}_bold.nii.gz'
		infile = infile[:-7]
		return infile

	def create_output_directory(self):
		output_dir = f'{self.output_dir}{self.subject_dir}/func'
		os.makedirs(output_dir, exist_ok=True)

	def identify_task_fmri_trimmed_file(self, task_type, run_number):
		outfile = f'{self.output_dir}{self.subject_dir}/func/{self.subject_dir}_task-{task_type}_run-0{run_number}_bold.nii.gz'
		outfile = outfile[:-7] + "_mod"
		return outfile

	def trim_task_fmri_volumes(self, number_volume=4):
		for i in range(1,4):
			infile = self.identify_task_fmri_input_file("EmotionRegulation", i)
			outfile = self.identify_task_fmri_trimmed_file("EmotionRegulation", i)
			os.system(f'fslroi {infile} {outfile} {number_volume} -1')
		
	def identify_resting_fmri_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		infile = infile[:-7]
		return infile

	def identify_resting_fmri_trimmed_file(self):
		outfile = f'{self.output_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		outfile = outfile[:-7] + "_mod"
		return outfile

	def trim_resting_fmri_volumes(self, number_volume=4):
		infile = self.identify_resting_fmri_input_file()
		outfile = self.identify_resting_fmri_trimmed_file()
		os.system(f'fslroi {infile} {outfile} {number_volume} -1')

	def process(self):
		self.create_output_directory()
		self.trim_task_fmri_volumes()
		self.trim_resting_fmri_volumes()		

class BiasCorrector:

	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	def identify_mask_output_file(self):
		mask_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected_brain_mask.nii.gz'
		return mask_outfile

	def create_tempory_directory(self):
		temporary_directory = f'{self.output_dir}{self.subject_dir}/anat/tmp/'
		os.makedirs(temporary_directory, exist_ok=True)

	def identify_temporary_directory(self):
		temporary_directory = f'{self.output_dir}{self.subject_dir}/anat/tmp/'
		return temporary_directory

	def identify_first_iteration_output_file(self):
		temporary_directory = self.identify_temporary_directory()
		first_iteration_outfile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter1.nii.gz'
		return first_iteration_outfile

	def identify_final_iteration_output_file(self):
		iteration_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected.nii.gz'
		return iteration_outfile

	def bias_correct(self):
		infile = self.identify_input_file()
		mask_outfile = self.identify_mask_output_file()
		first_iteration_outfile = self.identify_first_iteration_output_file()
		temporary_directory = self.identify_temporary_directory()
		os.system(f'N4BiasFieldCorrection -i {infile} -x {mask_outfile} -o {first_iteration_outfile}')

		for i in range(1, 6):
			iteration_infile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter{i}.nii.gz'
			i += 1 
			temporary_iteration_outfile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter{i}.nii.gz'
			os.system(f'N4BiasFieldCorrection -i {iteration_infile} -x {mask_outfile} -o {temporary_iteration_outfile}')

	def extract_final_iteration(self):
		temporary_directory = self.identify_temporary_directory()
		infile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter6.nii.gz'
		outfile = self.identify_final_iteration_output_file()
		shutil.copyfile(infile, outfile)

	def remove_temporary_directory(self):
		temporary_directory = self.identify_temporary_directory()
		shutil.rmtree(temporary_directory)


	def process(self):
		self.create_tempory_directory()
		self.bias_correct()
		self.extract_final_iteration()
		self.remove_temporary_directory()

class BSEBrainExtractor:

	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.study_name_all_caps = study_name.upper()

	# def identify_input_file(self):
	# 	infile = f'{self.input_dir}{study_name_all_caps}_Imaging/{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
	# 	return infile

	# def identify_bet_output_file(self):
	# 	outfile = f'{self.output_dir}{study_name_all_caps}_Imaging_Analysis/{self.subject_dir}/anat/{self.subject_dir}_T1w_brain.nii.gz'
	# 	return outfile

	# def create_output_dir(self):
	# 	outdir = f'{self.output_dir}{study_name_all_caps}_Imaging_Analysis/{self.subject_dir}/anat'
	# 	os.makedirs(outdir, exist_ok = True)


	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	def identify_bse_output_file(self):
		bse_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected_brain.nii.gz'
		return bse_outfile
 
	def identify_bse_mask_file(self):
		bse_binarized_mask_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected_brain_mask.nii.gz'
		return bse_binarized_mask_outfile

	def skull_strip_bse(self, infile, outfile):
		os.system(f'bse -i {infile} -o {outfile} --auto')#Skull-strip unfiltered T1w

	def binarize_bse_mask(self, infile, outfile):
		#Create binary brain mask using fslmaths
		os.system(f'fslmaths {infile} -thr 0.5 -bin {outfile}')#Create Binary mask

	def process(self):
		infile = self.identify_input_file()
		outfile = self.identify_bse_output_file()
		mask_outfile = self.identify_bse_mask_file()
		self.skull_strip_bse(infile, outfile)
		self.binarize_bse_mask(outfile, mask_outfile)


class BetBrainExtractor:

	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.study_name_all_caps = study_name.upper()

	# def identify_input_file(self):
	# 	infile = f'{self.input_dir}{study_name_all_caps}_Imaging/{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
	# 	return infile

	# def identify_bet_output_file(self):
	# 	outfile = f'{self.output_dir}{study_name_all_caps}_Imaging_Analysis/{self.subject_dir}/anat/{self.subject_dir}_T1w_brain.nii.gz'
	# 	return outfile

	# def create_output_dir(self):
	# 	outdir = f'{self.output_dir}{study_name_all_caps}_Imaging_Analysis/{self.subject_dir}/anat'
	# 	os.makedirs(outdir, exist_ok = True)


	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	def create_bet_output_directory(self):
		bet_output_dir = f'{self.output_dir}{self.subject_dir}/anat'
		os.makedirs(bet_output_dir, exist_ok = True)

	def identify_bet_output_file(self):
		bet_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_brain.nii.gz'
		return bet_outfile

	def skull_strip_bet(self, infile, outfile):
		os.system(f'bet {infile} {outfile} -m -f 0.3 -R')

	def process(self):
		infile = self.identify_input_file()
		self.create_bet_output_directory()
		outfile = self.identify_bet_output_file()
		self.skull_strip_bet(infile, outfile)



class OnsetCreator:

	# def __init__(self, subject_number, scan_eprime_by_run_dir="/study/midus3/processed_data/Temporary/Small/, study_name="midus3"):
	# 	self.scan_eprime_by_run_dir = scan_eprime_by_run_dir
	# 	self.study_name = study_name
	# 	self.subject_number = subject_number

	def __init__(self, study_name, subject_number):
		self.scan_eprime_by_run_dir = f"/study/{study_name}/processed_data/Temporary/Small/"
		self.study_name = study_name
		self.onset_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.subject_number = subject_number

	def identify_target_data(self, task_number):
		# Read in tsv in pandas dataframe
		IAPS_df = pd.read_csv(f'{self.scan_eprime_by_run_dir}sub-{self.subject_number}_IAPS-0{task_number}.tsv', sep='\t')
		faces_df = pd.read_csv(f'{self.scan_eprime_by_run_dir}sub-{self.subject_number}_faces-0{task_number}.tsv', sep='\t')

		return IAPS_df, faces_df

	def feature_engineer_data(self, df1, df2):
		# join the two data frames along rows
		merged_df = pd.concat([df1, df2])
		# Sort everything by onset column
		sorted_df = merged_df.sort_values('onset')
		# Remove unnecessary blocks column
		del sorted_df['Blocks']
		# Use loc and isnull to replace no responses to "n/a"
		sorted_df.loc[(sorted_df['Response_Time_Face'] == 0), 'Response_Time_Face'] = "n/a"
		sorted_df.loc[(sorted_df['response'].isnull()), 'response'] = "n/a"

		onset_df = sorted_df
		return onset_df 

	def create_onset_data(self, df, task_number):
		# Write the datafram into tsv format
		df.to_csv(f'{self.onset_dir}sub-{self.subject_number}/func/sub-{self.subject_number}_task-EmotionRegulation_run-0{task_number}_events.tsv', sep='\t', index=None) # index = None --> removes first column of the index which are created by to_csv Pandas dataframe

	def process(self):
		for i in range(1,4):
			IAPS_df, faces_df = self.identify_target_data(i)
			onset_df = self.feature_engineer_data(IAPS_df, faces_df)
			self.create_onset_data(onset_df, i)
			

class ScanEprimeDivider:

	# def __init__(self, subject_number, scan_eprime_dir="/study/midus3/processed_data/Temporary/Big/", 
	# 	scan_eprime_by_run_output_dir="/study/midus3/processed_data/Temporary/Small/", study_name="midus3"):
	# 	self.scan_eprime_dir = scan_eprime_raw_dir
	# 	self.scan_eprime_by_run_output_dir = scan_eprime_output_dir
	# 	self.study_name = study_name
	# 	self.subject_number = subject_number


	def __init__(self, study_name, subject_number):
		self.scan_eprime_dir = f"/study/{study_name}/processed_data/Temporary/Big/"
		self.scan_eprime_by_run_output_dir = f"/study/{study_name}/processed_data/Temporary/Small/"
		self.study_name = study_name
		self.subject_number = subject_number

	def set_column_headers(self):
		column_header_list = ['onset',
		'duration',
		'database',
		'Response_Time_Face',
		'stimulus',
		'correct',
		'valence',
		'valenceFollowing',
		'gender',
		'response',
		'face_correct_response',
		'sociality',
		'groups',
		'onset_trimmed',
		'Blocks']

		return column_header_list

		#Define column header names
	def create_csv_writer(self, outfile):
		column_header_list = self.set_column_headers()
		row_writer = csv.DictWriter(outfile, fieldnames=column_header_list, delimiter='\t', lineterminator='\n') 
		return row_writer

	def create_csv_reader(self, infile):
		row_reader= csv.DictReader(infile, delimiter='\t') # Dict Reader uses the header column names
		return row_reader

	def write_header(self, outfile):
		row_writer = self.create_csv_writer(outfile)
		row_writer.writeheader()

	def compute_and_write_IAPS_values(self, infile, outfile, task_number):
		row_reader = self.create_csv_reader(infile)
		row_writer = self.create_csv_writer(outfile)

		# Compute the values (onset, duration, response time etc.)
		for rows in row_reader:
			block_number = rows['Blocks']
			# Cast float to maintain the decimal points
			IAPS_onset_time = float(rows['IAPSPicture.OnsetTime'])/float(1000)
			IAPS_offset_time = float(rows['IAPSPicture.OffsetTime'])/float(1000)
			IAPS_TTL = float(rows['WaitForTTL.RTTime'])/float(1000)
			IAPS_onset_TTL_adjusted = IAPS_onset_time - IAPS_TTL
			IAPS_duration = float(IAPS_offset_time - IAPS_onset_time)
			IAPS_onset_time_trimmed = IAPS_onset_TTL_adjusted - 8#(float(rows['IAPSPicture.OnsetTime'])-8000)/float(1000)
			IAPS_number = rows['PictureFile'][:4]
			IAPS_valence = rows['valencecategory']
			IAPS_socilaity = rows['Sociality']
			ones = rows['Group']

			task_number = str(task_number)
			if task_number in block_number:

				row_writer.writerow({'onset': IAPS_onset_TTL_adjusted,
				'duration':IAPS_duration,
				'Response_Time_Face': "n/a",
				'database': "IAPS",
				'stimulus': IAPS_number,
				'correct': "n/a",
				'valence': IAPS_valence,
				'valenceFollowing': "n/a",
				'gender': "n/a",
				'response': "n/a",
				'face_correct_response': "n/a",
				'sociality': IAPS_socilaity,
				'groups': ones,
				'onset_trimmed': IAPS_onset_time_trimmed,
				'Blocks':rows['Blocks']})

	def compute_and_write_face_values(self, infile, outfile, task_number):
		row_reader = self.create_csv_reader(infile)
		row_writer = self.create_csv_writer(outfile)

		for rows in row_reader:
		
			block_number = rows['Blocks']
			# Cast float to maintain the decimal points
			face_onset_time = float(rows['Face.OnsetTime'])/float(1000) #CG
			face_offset_time = float(rows['Face.OffsetTime'])/float(1000)
			face_TTL = float(rows['WaitForTTL.RTTime'])/float(1000)
			face_onset_TTL_adjusted = face_onset_time - face_TTL
			face_onset_time_trimmed = face_onset_TTL_adjusted - 8 #(float(rows['Face.OnsetTime'])-8000)/float(1000)
			face_response_time = float(rows['Face.RT'])/float(1000) #CB
			face_duration = float(face_offset_time - face_onset_time)
			face_number = rows['FaceFile'][:2]
			face_gender = rows['Gender']
			face_response = rows['Face.RESP']
			face_correct_response = rows['Face.CRESP']
			face_correct = rows['Face.ACC']
			ones = rows['Group']
			IAPS_valence = rows['valencecategory']

			task_number = str(task_number)
			if task_number in block_number:
				if face_correct == '1':
					variable = 'Y'
				elif face_correct == '0':
					variable = 'N'
				row_writer.writerow({'onset': face_onset_TTL_adjusted,
				'duration': face_duration,
				'Response_Time_Face': face_response_time,
				'database': "faces",
				'stimulus': face_number,
				'correct': variable,
				'valence': "n/a",
				'valenceFollowing': IAPS_valence,
				'gender': face_gender,
				'response': face_response,
				'face_correct_response': face_correct_response,
				'sociality': "n/a",
				'groups': ones,
				'onset_trimmed': face_onset_time_trimmed,
				'Blocks':rows['Blocks']})

	def process_IAPS(self):
		for i in range(1,4):
			with open(f'{self.scan_eprime_dir}sub-{self.subject_number}_task-ER_events.tsv' , 'r') as infile, open(f'{self.scan_eprime_by_run_output_dir}sub-{self.subject_number}_IAPS-0{i}.tsv', 'w', newline="") as outfile:
				self.create_csv_reader(infile)
				self.create_csv_writer(outfile)
				self.write_header(outfile)
				self.compute_and_write_IAPS_values(infile, outfile, i)

	def process_faces(self):
		for i in range(1,4):
			with open(f'{self.scan_eprime_dir}sub-{self.subject_number}_task-ER_events.tsv' , 'r') as infile, open(f'{self.scan_eprime_by_run_output_dir}sub-{self.subject_number}_faces-0{i}.tsv', 'w', newline="") as outfile:
				self.create_csv_reader(infile)
				self.create_csv_writer(outfile)
				self.write_header(outfile)
				self.compute_and_write_face_values(infile, outfile, i)

	def process(self):
		self.process_IAPS()
		self.process_faces()


class ScanEprimeConverter:

	# def __init__(self, subject_number, scan_eprime_raw_dir="/study/midus3/raw-data/scan_eprime/data/", 
	# 	scan_eprime_output_dir="/study/midus3/processed_data/Temporary/Big/", study_name="midus3"):
	# 	self.scan_eprime_raw_dir = scan_eprime_raw_dir
	# 	self.scan_eprime_output_dir = scan_eprime_output_dir
	# 	self.study_name = study_name
	# 	self.subject_number = subject_number

	def __init__(self, study_name, subject_number):
		self.scan_eprime_raw_dir = f"/study/{study_name}/raw-data/scan_eprime/data/"
		self.scan_eprime_output_dir = f"/study/{study_name}/processed_data/Temporary/Big/"
		self.study_name = study_name
		self.subject_number = subject_number

	# def inspect_scan_eprime(self):
	# 	raw_scan_eprime_files = sorted(glob.glob(self.scan_eprime_raw_dir + "/midus3_order[1-2]_eyetracking_v[0-9][0-9]-[0-9][0-9][0-9]-[0-9]*.txt"))
	# 	for file in raw_scan_eprime_files:
	# if os.path.isfile(f'{self.scan_eprime_output_dir} sub-{self.subject_number}_task-ER_events.tsv'):
	# 		pass
	# 	elif not os.path.isfile(f'{self.scan_eprime_output_dir} sub-{self.subject_number}_task-ER_events.tsv'):


	def identify_target_raw_eprime(self):
		infile = f'{self.scan_eprime_raw_dir}{self.study_name}_order[1-2]_eyetracking_v0[0-2]-{self.subject_number}-{self.subject_number}.txt'
		return infile

	def set_outfile_name(self):
		outfile = f'{self.scan_eprime_output_dir}sub-{self.subject_number}_task-ER_events.tsv'

		return outfile

	def process(self):
		infile = self.identify_target_raw_eprime()
		outfile = self.set_outfile_name()
		os.system(f'eprime2tabfile {infile} > {outfile}')

class NiftiRotator:

	def __init__(self, study_name, subject_number):
		self.nifti_dir = f"/study/{study_name}/{study_name.upper()}_Imaging/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# def identify_target_subject(self):
	# 	niftis = glob.glob(self.nifti_dir + self.study_name + "_Imaging/sub-[0-9][0-9][0-9]")

	# 	for nifti in niftis:
	# 		subject_dir = nifti.split('/')[5]
	# 		subject_number = subject_dir[4:]

	def fix_task_fMRI_orientation(self, run_number, old_scan_name, new_scan_name):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{old_scan_name}_{run_number}_bold.nii.gz'
		if os.path.isfile(target_nifti):
			os.system(f'fslreorient2std {target_nifti} {self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{new_scan_name}_{run_number}_bold.nii.gz')
			os.remove(target_nifti)

	def rename_task_fMRI_json(self, run_number, old_scan_name, new_scan_name):
		target_json = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{old_scan_name}_{run_number}_bold.json'
		if os.path.isfile(target_json):
			os.rename(target_json, f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{new_scan_name}_{run_number}_bold.json')

	def fix_resting_fMRI_orientation(self):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		if os.path.isfile(target_nifti):
			os.system(f'fslreorient2std {target_nifti} {target_nifti}')

	def fix_dwi_orientation(self):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/dwi/{self.subject_dir}_dwi.nii.gz'
		if os.path.isfile(target_nifti):
			os.system(f'fslreorient2std {target_nifti} {target_nifti}')

	def process(self):
		self.fix_task_fMRI_orientation("run-01", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-01", "task-ER", "task-EmotionRegulation")
		self.fix_task_fMRI_orientation("run-02", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-02", "task-ER", "task-EmotionRegulation")
		self.fix_task_fMRI_orientation("run-03", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-03", "task-ER", "task-EmotionRegulation")
		self.fix_resting_fMRI_orientation()
		self.fix_dwi_orientation()



class NiftiConverter:

	def __init__(self, study_name, subject_number):
		
		self.nifti_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.raw_dir = f"/study/{study_name}/raw-data/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.subject_dicom_path = self.raw_dir + subject_number


	# def inspect_dicoms(self):
	# 	# raw DICOMs in a set
	# 	raw_dicoms = glob.glob(self.raw_dir + "[0-9][0-9][0-9]")
	# 	dicom_set = set([int(raw.split('/')[4][0:]) for raw in raw_dicoms]) # extracts INTEGERS w/o "0"
	# 	#dicoms.sort()
	# 	return dicom_set

	# def inspect_niftis(self):
	# 	# scan subjects that already have niftis
	# 	niftis = glob.glob(self.nifti_dir + self.study_name + "_Imaging/sub-[0-9][0-9][0-9]")
	# 	nifti_set = set([int(nifti.split('/')[5][4:7]) for nifti in niftis]) # in python the set is called hash table
	# 	#niftis.sort()
	# 	return nifti_set

	# def identify_target_subject(self):
	# 	dicom_set = self.inspect_dicoms()
	# 	nifti_set = self.inspect_niftis()

	# 	for dicom_number in dicom_set:
	# 		# cast string and zero-pad the subject numbers
	# 		dicom_number_string = str(dicom_number)
	# 		dicom_number_zero_padded = ["0" for i in range(3-len(dicom_number_string))] + [dicom_number_string]
	# 		subject_number = ''.join(dicom_number_zero_padded)
	# 		subject_dir = "sub-" + subject_number
	# 		subject_dicom_path = self.raw_dir + subject_number

	# 		if dicom_number in nifti_set:
	# 			pass

	# 		elif dicom_number not in nifti_set:
	# 			return subject_dir, subject_dicom_path, subject_number

	def extract_scan_info(self, scan_path):
		scan_name = scan_path.split('/')[6][6:]
		return scan_name

	def convert(self, scan_path, scan_type, subject_dir, scan_name):
		os.system(f'convert_dcm2niix {scan_path} {self.nifti_dir}{self.subject_dir}/{scan_type}/{self.subject_dir}_{scan_name}.nii.gz')

	def process(self):

		#subject_dir, subject_dicom_path, subject_number = self.identify_target_subject()

		scans = sorted(glob.glob(self.subject_dicom_path + "/dicoms/*/*.tgz"))

		for scan in scans:
			scan_name = self.extract_scan_info(scan)

			if scan_name == "ORIG_T1w": # should be 'ORIG_T1w' for scans with pure-filtered T1w's
				scan_type = "anat"
				new_scan_name = "T1w"
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-1_bold":
				scan_type = "func"
				new_scan_name = "task-EmotionRegulation_run-01_bold"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-2_bold":
				scan_type = "func"
				new_scan_name = "task-EmotionRegulation_run-02_bold"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-3_bold":
				scan_type = "func"
				new_scan_name = "task-EmotionRegulation_run-03_bold"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-rest_bold":
				scan_type = "func"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, scan_name)

			if scan_name == "dwi":
				scan_type = "dwi"
				os.makedirs(self.nifti_dir + self.subject_dir + "/dwi/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, scan_name)

			if scan_name == "WATER__fmap":
				scan_type = "fmap"
				os.makedirs(self.nifti_dir + self.subject_dir + "/fmap/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, "magnitude")

			if scan_name == "FieldMap__fmap":
				scan_type = "fmap"
				os.makedirs(self.nifti_dir + self.subject_dir + "/fmap/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, "fieldmap")

			if scan_name == "asl":
				scan_type = "asl"
				os.makedirs(self.nifti_dir + self.subject_dir + "/asl/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, scan_name)

