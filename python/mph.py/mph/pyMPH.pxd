

from libcpp.vector cimport vector
from libcpp.queue cimport priority_queue
from libcpp.pair cimport pair

cdef extern from "../../../mph/main.cpp":
	ctypedef struct GradedMatrix:
		vector[vector[pair[int, int]]] matrix
		vector[vector[int]] row_grades
		vector[vector[int]] column_grades

	ctypedef double input_t
		
	GradedMatrix presentation(vector[vector[input_t]] _points, vector[int] _metrics, vector[input_t] _max_metric_values, vector[int] _filters, int hom_dim) except +
	
	GradedMatrix presentation_dm(vector[vector[vector[input_t]]] distance_matrices, vector[input_t] max_metric_values, vector[vector[input_t]] filters, int hom_dim) except +

	GradedMatrix presentation_FIrep(vector[vector[int]] high_matrix, vector[vector[int]] column_grades_h, vector[vector[int]] low_matrix, vector[vector[int]] column_grades_l) except +

	pair[GradedMatrix, GradedMatrix] groebner_bases(vector[vector[int]] matrix, vector[vector[int]] row_grades, vector[vector[int]] column_grades) except +

