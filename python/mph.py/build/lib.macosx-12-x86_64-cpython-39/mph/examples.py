import numpy as np
from typing import List, Tuple

from mph import GradedMatrix, groebner_bases, presentation_FIrep

def choose_graded_subbasis(matrix: List[List[Tuple[int, int]]], column_grades: List[List[int]], row_grades: List[List[int]]):
	dense_matrix = np.zeros(shape=(len(row_grades), len(column_grades)), dtype=np.int32)
	for i in range(len(matrix)):
		for entry in matrix[i]:
			dense_matrix[entry[1], i] = entry[0]
	transform = np.random.choice([0, 1], size=(dense_matrix.shape[1], np.random.randint(1, dense_matrix.shape[1]-1) if dense_matrix.shape[1]>2 else 1))
	column_grades_ret = []
	for i in range(transform.shape[1]):
		grade = [0]*len(column_grades[0])
		for j in range(transform.shape[0]):
			if transform[j, i] != 0:
				grade = [max(grade[k], column_grades[j][k]) for k in range(len(grade))]
		column_grades_ret.append(grade)
	return np.matmul(dense_matrix, transform) % 2, column_grades_ret


def generate_rivet_output(high_matrix, column_grades_h, low_matrix, column_grades_l):
	with open("input_file", "w") as f:
		out = "firep\nparameter 1\nparameter 2\n" + str(high_matrix.shape[1]) +" " + str(low_matrix.shape[1]) + " " + str(low_matrix.shape[0])+"\n"
		for column_index in range(high_matrix.shape[1]):
			out += str(column_grades_h[column_index][0]) +" "+str(column_grades_h[column_index][1])+" ; "
			for row_index in range(high_matrix.shape[0]):
				if high_matrix[row_index, column_index] != 0:
					out += str(row_index)+" "
			out += "\n"
		for column_index in range(low_matrix.shape[1]):
			out += str(column_grades_l[column_index][0]) +" "+str(column_grades_l[column_index][1])+" ; "
			for row_index in range(low_matrix.shape[0]):
				if low_matrix[row_index, column_index] != 0:
					out += str(row_index)+" "
			out += "\n"
		f.write(out)



def rivet_benchmark(m: int, n: int, n_parameters: int=2, density: float=0.5, grade_range:int=1000, log_level:str="silent"):
	# Create a random matrix as the low_matrix in the FIrep
	low_matrix = np.random.choice([0, 1], size=(m, n), p=[density, 1-density])
	# Choose random grades for the columns
	column_grades_l = [ np.random.choice(grade_range, size=(n_parameters,)) for _ in range(n) ]
	# Compute the kernel of this graded matrix
	output = groebner_bases(low_matrix, column_grades_l, log_level=log_level)
	if len(output[1].column_grades) == 0:
		return 0
	# Choose a random subbasis of the kernel as the high_matrix in the FIrep
	high_matrix, column_grades_h = choose_graded_subbasis(output[1].matrix, output[1].column_grades, output[1].row_grades)
	
	generate_rivet_output(high_matrix, column_grades_h, low_matrix, column_grades_l)

	ret = presentation_FIrep(high_matrix, column_grades_h, low_matrix, column_grades_l, log_level=log_level)

def random_FIrep_presentation(m: int, n: int, n_parameters: int=3, density: float=0.5, grade_range:int=100):
	# Create a random matrix as the low_matrix in the FIrep
	low_matrix = np.random.choice([0, 1], size=(m, n), p=[density, 1-density])
	# Choose random grades for the columns
	column_grades_l = [ np.random.choice(grade_range, size=(n_parameters,)) for _ in range(n) ]
	# Compute the kernel of this graded matrix
	output = groebner_bases(low_matrix, column_grades_l)
	if len(output[1].column_grades) == 0:
		return 0
	# Choose a random subbasis of the kernel as the high_matrix in the FIrep
	high_matrix, column_grades_h = choose_graded_subbasis(output[1].matrix, output[1].column_grades, output[1].row_grades)
	return presentation_FIrep(high_matrix, column_grades_h, low_matrix, column_grades_l)

def random_map_gbs(m:int, n:int, n_parameters: int=3, density:float=0.5, grade_range:int=100):
	random_matrix = np.random.choice([0, 1], size=(m, n), p=[density, 1-density])
	# Choose random grades for the columns
	column_grades = [ np.random.choice(grade_range, size=(n_parameters,)) for _ in range(n) ]
	# Compute Groebner bases for the image and kernel of this graded matrix
	image_gb, kernel_gb = groebner_bases(random_matrix, column_grades)
	sparse_random_matrix = [[] for _ in range(random_matrix.shape[1])]
	non_zero_entries = np.nonzero(random_matrix)
	for entry in list(zip(non_zero_entries[1], non_zero_entries[0])):
		sparse_random_matrix[entry[0]].append((random_matrix[entry[1], entry[0]], entry[1]))
	return GradedMatrix(sparse_random_matrix, column_grades, row_grades = [ [ 0 for _ in range(len(column_grades[0])) ] for _ in range(random_matrix.shape[0])]), image_gb, kernel_gb

