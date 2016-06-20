import numpy
cimport numpy

cimport cython

from libc.stdio cimport fscanf, fopen

@cython.boundscheck(False)
def build_label_array_from_bed(labels_bed):
    fname = labels_bed.fn
    cdef int i, start, stop
    i = 0
    ptr_fr = fopen(fname, "rb")
    num_lines = labels_bed.count()
    cdef int[::1] ptr_labels = numpy.empty((num_lines,), dtype='int32')
    for i in xrange(num_lines):
        fscanf(ptr_fr,"%*s\t%i\t%i\t%i\n", &start, &stop, &ptr_labels[i])
        #print start, stop, ptr_labels[i]
        i += 1
    return numpy.asarray(ptr_labels)
