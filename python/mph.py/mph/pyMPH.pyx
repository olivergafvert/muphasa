# distutils: language = c++

cimport numpy as np
cimport mph.pyMPH as pyMPH
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def doPresentation(np.ndarray[double,ndim=2,mode="c"] points not None, np.ndarray[int,ndim=1,mode="c"] metrics, np.ndarray[double,ndim=1,mode="c"] max_metric_values, np.ndarray[int,ndim=1,mode="c"] filters, int hom_dim):

	res = pyMPH.presentation(points, metrics, max_metric_values, filters, hom_dim)

	return res

@cython.boundscheck(False)
@cython.wraparound(False)
def doPresentationDm(np.ndarray[double, ndim=3, mode="c"] distance_matrices not None, np.ndarray[double, ndim=1, mode="c"] max_metric_values, np.ndarray[double, ndim=2, mode="c"] filters, int hom_dim):

	res = pyMPH.presentation_dm(distance_matrices, max_metric_values, filters, hom_dim)

	return res

@cython.boundscheck(False)
@cython.wraparound(False)
def doPresentationFIrep(np.ndarray[int,ndim=2,mode="c"] high_matrix not None, np.ndarray[int,ndim=2,mode="c"] column_grades_h, np.ndarray[int,ndim=2,mode="c"] low_matrix not None, np.ndarray[int,ndim=2,mode="c"] column_grades_l):

	res = pyMPH.presentation_FIrep(high_matrix, column_grades_h, low_matrix, column_grades_l)

	return res

@cython.boundscheck(False)
@cython.wraparound(False)
def computeGroebnerBases(np.ndarray[int,ndim=2,mode="c"] matrix not None, np.ndarray[int,ndim=2,mode="c"] row_grades, np.ndarray[int,ndim=2,mode="c"] column_grades):

	res = pyMPH.groebner_bases(matrix, row_grades, column_grades)

	return res

