#!/usr/bin/env python

__author__ = "Seung Beum (David) Lee"
__email__ = "sleedave90@gmail.com"
__credits__ = ["Seung Beum (David) Lee", "Micheal Kelly", "Jeanette Mumford", "Nate Vack"]

import os
import glob
import csv
import pandas as pd
import shutil
import os.path
import json

################ READ this before editing brainiac ################
# Each class serves distinct purpose 
# Each class has "process" function runs all subsequent functions
# Each class is imported in a wraper script 
# Each class sets directory paths
# If the program is used for datasets other than MIDUS3 / EMOWRAP, the paths need to be re-coded. 
# All classes and functions are named after their purpose

################ Global Variables ################ 
# Edit these variables if needed
number_volumes_trimmed = 4 # How many volumes (TRs) of fMRI data do you want to trim to account for scanner warm up?
number_biascorrection_iteration = 6 # How many iterations of N4biascorrection do you want to perform?
fractional_intensity_threshold_for_brain_extraction = 0.3 # Value between (0->1); default=0.5; smaller values give larger brain outline estimates
task_fmri_framewise_displacement_threshold = 0.5 # unit = mm
resting_fmri_framewise_displacement_threshold = 0.2 # unit = mm, may change to 0.3 if too strict

# Create BIDS compatible JSON file that is required to pass BIDS vaidator
# Passing the validator is required to use any of BIDS compatible apps
class JsonCreator:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# Set output json file name
	def identify_json(self):
		json = f'{self.input_dir}dataset_description.json'
		return json

	# Check existence of a json file
	def check_json(self, json_dict):
		json = self.identify_json()
		
		try:
			os.system(f"more {json}")
			print ("=====Dataset json file exists, check content above=====")
		except IOError:
			with open (json, 'w') as file:
				file.write(json.dumps(json_dict))
			print ("=====Dataset json file does not exist, creating one=====")

	# Create a json file if it doesn't exist
	def create_json(self):
		json_dict = {}
		key_list = ["Name", "BIDSVersion", "License", "Authors", "Acknowledgements", "HowToAcknowledge", "Funding", "ReferencesAndLinks", "DatasetDOI"]
		for key in key_list:
			value = input(f"Enter Value for {key}: ")
			json_dict[key] = value

		print ("")
		print (json_dict)
		return (json_dict)

	# process all
	def process(self):
		self.check_json(self)

# Assess motion in all fMRI data
class MotionEvaluator:

	# Set paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# Set output directory name
	def identify_qa_directory(self):
		qa_dir = f'{self.output_dir}{self.subject_dir}/func/QA'
		return qa_dir

	# Create output directory in derivative data folder
	def create_qa_directory(self):
		qa_dir = self.identify_qa_directory()
		os.makedirs(qa_dir, exist_ok=True)

	# Set output html file name that will contain motion assessment results
	def identify_html_file(self):
		qa_dir = self.identify_qa_directory()
		html_file = f'{qa_dir}/motion.html'
		return html_file

	# Perform motion assessment in all runs of task fMRI data
	def evaluate_task_fmri_motion(self):
		qa_dir = self.identify_qa_directory()
		html_file = self.identify_html_file()
		volume_trimmer = VolumeTrimmer(self.study_name, self.subject_number)
		for i in range(1,4):
			infile = volume_trimmer.identify_task_fmri_trimmed_file("EmotionRegulation", i)
			confound_outfile = f'{qa_dir}/0{i}_confound.txt'
			plot_outfile = f'{qa_dir}/fd_plot_0{i}.png'
			outlier_outfile = f'{qa_dir}/outlier_output_0{i}.txt'
			print (f"=====assesing motion for task-fMRI run {i}=====")
			os.system(f'fsl_motion_outliers -i {infile} -o {confound_outfile} --fd --thresh={task_fmri_framewise_displacement_threshold} -p {plot_outfile} -v > {outlier_outfile}')
			if not os.path.isfile(confound_outfile):
				os.mknod(confound_outfile)
			os.system(f'cat {outlier_outfile} >> {html_file}')
			os.system(f"echo '<br><br> FD plot {self.subject_number} 0{i} <br><IMG BORDER=0 SRC={plot_outfile} WIDTH=100%></BODY></HTML> <p>==========================================================<p>' >> {html_file}")

	# Perform motion assessment in resting fMRI data
	def evaluate_resting_fmri_motion(self):
		qa_dir = self.identify_qa_directory()
		html_file = self.identify_html_file()
		volume_trimmer = VolumeTrimmer(self.study_name, self.subject_number)
		infile = volume_trimmer.identify_resting_fmri_trimmed_file()
		confound_outfile = f'{qa_dir}/04_confound.txt'
		plot_outfile = f'{qa_dir}/fd_plot_04.png'
		outlier_outfile = f'{qa_dir}/outlier_output_04.txt'
		print (f"=====assesing motion for resting-fMRI=====")
		os.system(f'fsl_motion_outliers -i {infile} -o {confound_outfile} --fd --thresh={resting_fmri_framewise_displacement_threshold} -p {plot_outfile} -v > {outlier_outfile}')
		if not os.path.isfile(confound_outfile):
				os.mknod(confound_outfile)
		os.system(f'cat {outlier_outfile} >> {html_file}')
		os.system(f"echo '<br><br> FD plot {self.subject_number} 04 <br><IMG BORDER=0 SRC={plot_outfile} WIDTH=100%></BODY></HTML> <p>==========================================================<p>' >> {html_file}")

	# Process all
	def process(self):
		self.create_qa_directory()
		self.evaluate_task_fmri_motion()
		self.evaluate_resting_fmri_motion()

# Trim first 4 volumes of fMRI data to account for scanner warm-up
class VolumeTrimmer:

	# Set paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# Identify input task fMRI data
	def identify_task_fmri_input_file(self, task_type, run_number):
		infile = f'{self.input_dir}{self.subject_dir}/func/{self.subject_dir}_task-{task_type}_run-0{run_number}_bold.nii.gz'
		infile = infile[:-7]
		return infile

	# Create ouput directory in derivative data folder
	def create_output_directory(self):
		output_dir = f'{self.output_dir}{self.subject_dir}/func'
		os.makedirs(output_dir, exist_ok=True)

	# Set output name for trimmed task fMRI data
	def identify_task_fmri_trimmed_file(self, task_type, run_number):
		outfile = f'{self.output_dir}{self.subject_dir}/func/{self.subject_dir}_task-{task_type}_run-0{run_number}_bold.nii.gz'
		outfile = outfile[:-7] + "_mod"
		return outfile

	# Trim first 4 volumes of all runs of task fMRI data
	def trim_task_fmri_volumes(self, number_volume=number_volumes_trimmed):
		for i in range(1,4):
			infile = self.identify_task_fmri_input_file("EmotionRegulation", i)
			outfile = self.identify_task_fmri_trimmed_file("EmotionRegulation", i)
			print (f"=====trimming first {str(number_volume)} volmues of run {i}=====")
			os.system(f'fslroi {infile} {outfile} {number_volume} -1')
		
	# Identify input resting fMRI data
	def identify_resting_fmri_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		infile = infile[:-7]
		return infile

	# Set output name for trimmed resting fMRI data
	def identify_resting_fmri_trimmed_file(self):
		outfile = f'{self.output_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		outfile = outfile[:-7] + "_mod"
		return outfile

	# Trim first 4volumes of resting fMRI data
	def trim_resting_fmri_volumes(self, number_volume=4):
		infile = self.identify_resting_fmri_input_file()
		outfile = self.identify_resting_fmri_trimmed_file()
		os.system(f'fslroi {infile} {outfile} {number_volume} -1')
		print (f"=====trimming first {str(number_volume)} volmues of resting-fMRI=====")

	# Process all
	def process(self):
		self.create_output_directory()
		self.trim_task_fmri_volumes()
		self.trim_resting_fmri_volumes()		

# Biascorrect raw T1w images using ANTs N4biascorrection
class BiasCorrector:

	# Set paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# Identify raw T1w input iamge
	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	# Set output file name
	def identify_mask_output_file(self, mask_input):
		mask_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_{mask_input}'#T1w_brain_mask.nii.gz
		return mask_outfile

	# Create temporay directory that will contatin every iterations of biascorrected images
	def create_tempory_directory(self):
		temporary_directory = self.identify_temporary_directory()
		os.makedirs(temporary_directory, exist_ok=True)

	# Set temporary directory name
	def identify_temporary_directory(self):
		temporary_directory = f'{self.output_dir}{self.subject_dir}/anat/tmp/'
		return temporary_directory

	# Set first-iteration output image name
	def identify_first_iteration_output_file(self):
		temporary_directory = self.identify_temporary_directory()
		first_iteration_outfile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter1.nii.gz'
		return first_iteration_outfile

	# Set input image name
	def identify_first_final_iteration_output_file(self):
		iteration_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected.nii.gz'
		return iteration_outfile

	# Set ouput image name
	def identify_second_final_iteration_output_file(self):
		iteration_second_outfile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected.nii.gz'
		return iteration_second_outfile

	# Perform bias correction
	def bias_correct(self, mask_input):
		infile = self.identify_input_file()
		mask_outfile = self.identify_mask_output_file(mask_input)
		first_iteration_outfile = self.identify_first_iteration_output_file()
		temporary_directory = self.identify_temporary_directory()
		os.system(f'N4BiasFieldCorrection -i {infile} -x {mask_outfile} -o {first_iteration_outfile}')

		for i in range(1, number_biascorrection_iteration):
			iteration_infile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter{i}.nii.gz'
			i += 1 
			temporary_iteration_outfile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter{i}.nii.gz'
			os.system(f'N4BiasFieldCorrection -i {iteration_infile} -x {mask_outfile} -o {temporary_iteration_outfile}')

	# Identify final product
	def extract_final_iteration(self):
		temporary_directory = self.identify_temporary_directory()
		infile = f'{temporary_directory}{self.subject_dir}_T1w_N4Corrected_iter6.nii.gz'
		first_outfile = self.identify_first_final_iteration_output_file()
		second_outfile = self.identify_second_final_iteration_output_file()
		shutil.copyfile(infile, first_outfile)
		shutil.copyfile(infile, second_outfile)

	# Remove temporary directory with non-usable iterative images
	def remove_temporary_directory(self):
		temporary_directory = self.identify_temporary_directory()
		shutil.rmtree(temporary_directory)

	# Process all
	def process(self, mask_input = 'T1w_N4Corrected_brain_mask.nii.gz'):
		self.create_tempory_directory()
		self.bias_correct(mask_input)
		self.extract_final_iteration()
		self.remove_temporary_directory()

# Skull-strip using BSE
class BSEBrainExtractor:

	# Set paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.study_name_all_caps = study_name.upper()

	# identify T1w image as input
	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	# set ouput name
	def identify_bse_output_file(self):
		bse_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected_brain.nii.gz'
		return bse_outfile
 
 	# set ouput binary mask file name
	def identify_bse_mask_file(self):
		bse_binarized_mask_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_N4Corrected_brain_mask.nii.gz'
		return bse_binarized_mask_outfile

	# Perform skull stripping (BSE had to be installed through BIT help)
	def skull_strip_bse(self, infile, outfile):
		os.system(f'bse -i {infile} -o {outfile} --auto')#Skull-strip unfiltered T1w

	# binarize brain mask file
	def binarize_bse_mask(self, infile, outfile):
		#Create binary brain mask using fslmaths
		os.system(f'fslmaths {infile} -thr 0.5 -bin {outfile}')#Create Binary mask

	# process all
	def process(self):
		infile = self.identify_input_file()
		outfile = self.identify_bse_output_file()
		mask_outfile = self.identify_bse_mask_file()
		self.skull_strip_bse(infile, outfile)
		self.binarize_bse_mask(outfile, mask_outfile)

# Skull-strip using fsl bet
class BetBrainExtractor:

	# Set paths
	def __init__(self, study_name, subject_number):
		self.input_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.output_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging_Analysis/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.study_name_all_caps = study_name.upper()

	# identify T1w image as input
	def identify_input_file(self):
		infile = f'{self.input_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		return infile

	# create "anat" directory in a derivative folder
	def create_bet_output_directory(self):
		bet_output_dir = f'{self.output_dir}{self.subject_dir}/anat'
		os.makedirs(bet_output_dir, exist_ok = True)

	# set ouput name
	def identify_bet_output_file(self):
		bet_outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_brain.nii.gz'
		return bet_outfile

	def skull_strip_bet(self, infile, outfile):
		os.system(f'bet {infile} {outfile} -m -f {fractional_intensity_threshold_for_brain_extraction} -R')

	# Copy T1w image in to derivative folder
	def copy_original_T1w(self):
		infile = self.identify_input_file()
		outfile = f'{self.output_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		shutil.copyfile(infile, outfile)

	# Process all
	def process(self):
		infile = self.identify_input_file()
		self.create_bet_output_directory()
		outfile = self.identify_bet_output_file()
		self.skull_strip_bet(infile, outfile)
		self.copy_original_T1w()


# Create FSL compatible 3-column timing files for future task fMRI analysis
class OnsetCreator:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.scan_eprime_by_run_dir = f"/study/{study_name}/processed_data/Temporary/Small/"
		self.study_name = study_name
		self.onset_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.subject_number = subject_number

	# Identify E-prime files for each run
	def identify_target_data(self, task_number):
		# Read in tsv in pandas dataframe
		IAPS_df = pd.read_csv(f'{self.scan_eprime_by_run_dir}sub-{self.subject_number}_IAPS-0{task_number}.tsv', sep='\t')
		faces_df = pd.read_csv(f'{self.scan_eprime_by_run_dir}sub-{self.subject_number}_faces-0{task_number}.tsv', sep='\t')

		return IAPS_df, faces_df

	# Merge IAPS/neutral face EPrime files 
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

	# Ouput to .csv
	def create_onset_data(self, df, task_number):
		# Write the datafram into tsv format
		df.to_csv(f'{self.onset_dir}sub-{self.subject_number}/func/sub-{self.subject_number}_task-EmotionRegulation_run-0{task_number}_events.tsv', sep='\t', index=None) # index = None --> removes first column of the index which are created by to_csv Pandas dataframe

	# Process all
	def process(self):
		for i in range(1,4):
			IAPS_df, faces_df = self.identify_target_data(i)
			onset_df = self.feature_engineer_data(IAPS_df, faces_df)
			self.create_onset_data(onset_df, i)
			
# Extract columns of interest from giant Eprime output - columns that need be QA'd and used for future analyses
class ScanEprimeDivider:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.scan_eprime_dir = f"/study/{study_name}/processed_data/Temporary/Big/"
		self.scan_eprime_by_run_output_dir = f"/study/{study_name}/processed_data/Temporary/Small/"
		self.study_name = study_name
		self.subject_number = subject_number

	# Set column header names
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

	# Create a CSV writer with a list of column names
	def create_csv_writer(self, outfile):
		column_header_list = self.set_column_headers()
		row_writer = csv.DictWriter(outfile, fieldnames=column_header_list, delimiter='\t', lineterminator='\n') 
		return row_writer

	# Create a CSV reader 
	def create_csv_reader(self, infile):
		row_reader= csv.DictReader(infile, delimiter='\t') # Dict Reader uses the header column names
		return row_reader

	# Create a CSV header writer 
	def write_header(self, outfile):
		row_writer = self.create_csv_writer(outfile)
		row_writer.writeheader()

	# Initialize CSV reader and writer, then extract columns of interest from the original E-Prime file (IAPS trials)
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

	# Initialize CSV reader and writer, then extract columns of interest from the original E-Prime file (neutral face trials)
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

	# Process IAPS trials
	def process_IAPS(self):
		for i in range(1,4):
			with open(f'{self.scan_eprime_dir}sub-{self.subject_number}_task-ER_events.tsv' , 'r') as infile, open(f'{self.scan_eprime_by_run_output_dir}sub-{self.subject_number}_IAPS-0{i}.tsv', 'w', newline="") as outfile:
				self.create_csv_reader(infile)
				self.create_csv_writer(outfile)
				self.write_header(outfile)
				self.compute_and_write_IAPS_values(infile, outfile, i)

	# Process neutral face trials
	def process_faces(self):
		for i in range(1,4):
			with open(f'{self.scan_eprime_dir}sub-{self.subject_number}_task-ER_events.tsv' , 'r') as infile, open(f'{self.scan_eprime_by_run_output_dir}sub-{self.subject_number}_faces-0{i}.tsv', 'w', newline="") as outfile:
				self.create_csv_reader(infile)
				self.create_csv_writer(outfile)
				self.write_header(outfile)
				self.compute_and_write_face_values(infile, outfile, i)

	# Process all column extraction 
	def process(self):
		self.process_IAPS()
		self.process_faces()


class ScanEprimeConverter:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.scan_eprime_raw_dir = f"/study/{study_name}/raw-data/scan_eprime/data/"
		self.scan_eprime_output_dir = f"/study/{study_name}/processed_data/Temporary/Big/"
		self.study_name = study_name
		self.subject_number = subject_number

	# Identify raw eprime data
	def identify_target_raw_eprime(self):
		infile = f'{self.scan_eprime_raw_dir}{self.study_name}_order[1-2]_eyetracking_v0[0-2]-{self.subject_number}-{self.subject_number}.txt'
		return infile

	# Set output name
	def set_outfile_name(self):
		outfile = f'{self.scan_eprime_output_dir}sub-{self.subject_number}_task-ER_events.tsv'

		return outfile

	# Convert .txt raw eprime into .tsv using the latest conversion tool on the server
	def process(self):
		infile = self.identify_target_raw_eprime()
		outfile = self.set_outfile_name()
		os.system(f'eprime2tabfile {infile} > {outfile}')

# Re-orient T1w and fMRI data 
class NiftiRotator:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.nifti_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number

	# Re-orient T1w
	def fix_T1w_orientation(self):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w.nii.gz'
		if os.path.isfile(target_nifti):
			print (f"=====fixing orientations for {self.subject_number} T1w=====")
			os.system(f'fslreorient2std {target_nifti} {self.nifti_dir}{self.subject_dir}/anat/{self.subject_dir}_T1w_reoriented.nii.gz')

	# Re-orient task fMRI 
	# Then remove the old file
	def fix_task_fMRI_orientation(self, run_number, old_scan_name, new_scan_name):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{old_scan_name}_{run_number}_bold.nii.gz'
		#print (target_nifti)
		if os.path.isfile(target_nifti):
			print (f"=====fixing orientations for {self.subject_number} task-fMRI {run_number}=====")
			os.system(f'fslreorient2std {target_nifti} {self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{new_scan_name}_{run_number}_bold.nii.gz')
			os.remove(target_nifti)

	# Re-name task fMRI json file
	def rename_task_fMRI_json(self, run_number, old_scan_name, new_scan_name):
		target_json = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{old_scan_name}_{run_number}_bold.json'
		if os.path.isfile(target_json):
			os.rename(target_json, f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_{new_scan_name}_{run_number}_bold.json')

	# Re-oreint resting fMRI
	# Then remove the old file
	def fix_resting_fMRI_orientation(self):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/func/{self.subject_dir}_task-rest_bold.nii.gz'
		if os.path.isfile(target_nifti):
			print (f"=====fixing orientations for {self.subject_number} resting-fMRI=====")
			os.system(f'fslreorient2std {target_nifti} {target_nifti}')

	# Re-orient DWI data
	def fix_dwi_orientation(self):
		target_nifti = f'{self.nifti_dir}{self.subject_dir}/dwi/{self.subject_dir}_dwi.nii.gz'
		if os.path.isfile(target_nifti):
			print (f"=====fixing orientations for {self.subject_number} DWI=====")
			os.system(f'fslreorient2std {target_nifti} {target_nifti}')

	# Perform all re-orientation
	def process(self):
		self.fix_T1w_orientation()
		self.fix_task_fMRI_orientation("run-01", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-01", "task-ER", "task-EmotionRegulation")
		self.fix_task_fMRI_orientation("run-02", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-02", "task-ER", "task-EmotionRegulation")
		self.fix_task_fMRI_orientation("run-03", "task-ER", "task-EmotionRegulation")
		self.rename_task_fMRI_json("run-03", "task-ER", "task-EmotionRegulation")
		self.fix_resting_fMRI_orientation()
		self.fix_dwi_orientation()


# Convert DICOMs to NIFTIs
class DicomConverter:

	# Set Paths
	def __init__(self, study_name, subject_number):
		self.nifti_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.raw_dir = f"/study/{study_name}/raw-data/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.subject_dicom_path = self.raw_dir + subject_number

	# From current raw-data directory structure, extract the full scan name (scan name was explicitly determined befroe during piloting - also can be changed anytime from scanner protocol)
	def extract_scan_info(self, scan_path):
		scan_name = scan_path.split('/')[6][6:]
		return scan_name

	# Use the latest conversion tool from the server (takes zipped .tgz dicoms as input and creates NIFTI)
	def convert(self, scan_path, scan_type, subject_dir, scan_name):
		os.system(f'convert_dcm2niix {scan_path} {self.nifti_dir}{self.subject_dir}/{scan_type}/{self.subject_dir}_{scan_name}.nii.gz')

	# Hard-coded each scan paramter
	# Add in any additional imaging modality in the future
	def process(self):

		#subject_dir, subject_dicom_path, subject_number = self.identify_target_subject()

		scans = sorted(glob.glob(self.subject_dicom_path + "/dicoms/*/*.tgz"))

		for scan in scans:
			scan_name = self.extract_scan_info(scan)

			if scan_name == "ORIG_T1w": # should be 'ORIG_T1w' for scans with pure-filtered T1w's
				scan_type = "anat" # this is the directory name 
				new_scan_name = "T1w" # this is the file name
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "CUBE_T2_FLAIR": # This is the pure filtered flair Image (used for testing purposes, may have to remove if proven useless)
				scan_type = "anat"
				new_scan_name = "T2_flair_pure_filtered"
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "CUBE_T2": # This is the pure filtered cube Image (used for testing purposes, may have to remove if proven useless)
				scan_type = "anat" 
				new_scan_name = "T2_cube_pure_filtered"  
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "ORIG_CUBE_T2_FLAIR": 
				scan_type = "anat" 
				new_scan_name = "T2_flair" 
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "ORIG_CUBE_T2":
				scan_type = "anat"
				new_scan_name = "T2_cube"
				os.makedirs(self.nifti_dir + self.subject_dir + "/anat/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-1_bold":
				scan_type = "func"
				new_scan_name = "task-ER_run-01_bold"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-2_bold":
				scan_type = "func"
				new_scan_name = "task-ER_run-02_bold"
				os.makedirs(self.nifti_dir + self.subject_dir + "/func/", exist_ok=True)
				self.convert(scan, scan_type, self.subject_dir, new_scan_name)

			if scan_name == "task-ER_run-3_bold":
				scan_type = "func"
				new_scan_name = "task-ER_run-03_bold"
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

# Checks for duplicates to prevent overwriting
class ProcessChecker:

	# Set Paths
	def __init__(self, study_name, subject_number):
		
		self.nifti_dir = f"/study/{study_name}/processed_data/{study_name.upper()}_Imaging/"
		self.raw_dir = f"/study/{study_name}/raw-data/"
		self.study_name = study_name
		self.subject_number = subject_number
		self.subject_dir = "sub-" + subject_number
		self.subject_dicom_path = self.raw_dir + subject_number
		self.subject_nifti_path = self.nifti_dir + self.subject_dir

	# Check existence of raw-data
	def check_raw_data(self):
		if os.path.isdir(self.subject_dicom_path):
			return True

		else:
			return False

	# Check existence of processed data
	def check_processed_data(self):
		if os.path.isdir(self.subject_nifti_path):
			return True
		else:
			return False

	# Evaluate every scenario and outputs recommended course of action
	def process(self):

		check_raw_boolean = self.check_raw_data()
		check_processed_data_boolean = self.check_processed_data()

		if check_raw_boolean == True and check_processed_data_boolean == False:
			print (f"=====Niftis for {self.subject_number} do not exist, INITIATE processing=====")

		elif check_raw_boolean == True and check_processed_data_boolean == True:
			print (f"=====Raw Dicoms and Processed Niftis for {self.subject_number} already exist, ABORT processing=====")
			exit()

		elif check_raw_boolean == False and check_processed_data_boolean == True:
			print (f"=====Raw Dicoms for {self.subject_number} do not exist, but porcessed Niftis do, RED FLAG, check data=====")
			exit()

		elif check_raw_boolean == False and check_processed_data_boolean == False:
			print (f"=====Raw Dicoms and Processed Niftis for {self.subject_number} do not exist, check study name and subject number=====")
			exit()
