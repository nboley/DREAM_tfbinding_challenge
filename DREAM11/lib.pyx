import numpy
cimport numpy

cimport cython

from libc.stdio cimport fscanf, fopen, fclose, fgetc, EOF

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

@cython.boundscheck(False)
def cmp_labels(truth_bed_fname, submitted_bed_fname):
    """Returns true if the labels are identical in truth_bed_fname and submitted_bed_fname

    """
    cdef int truth_scanf_rv, submitted_scanf_rv
    cdef int truth_start, truth_stop
    cdef int submitted_start, submitted_stop
    cdef int truth_num_lines, submitted_num_lines
    cdef char[256] submitted_contig
    cdef char[256] truth_contig
    
    i = 0
    truth_fp = fopen(truth_bed_fname, "rb")
    submitted_fp = fopen(submitted_bed_fname, "rb")
    ptr_fr = fopen(fname, "rb")
    while True:
        fscanf(ptr_fr,"%s\t%i\t%i\n", &start, &stop, &ptr_labels[i])
        #print start, stop, ptr_labels[i]
        i += 1
    return numpy.asarray(ptr_labels)

@cython.boundscheck(False)
def load_labels(fname, max_num_lines):
    cdef int num_factors, i, j, start, stop, num_lines
    cdef char ch
    
    # load the header, and find the number of factors
    with open(fname) as fp:
        header_data = next(iter(fp)).split()
        if header_data[:3] != ['chr', 'start', 'stop']:
            raise ValueError, "Unexpected header"
        num_factors = len(header_data) - 3
    
    # find the number of lines
    ptr_fr = fopen(fp.name, "rb")
    num_lines = 1
    while True:
        ch = fgetc(ptr_fr)
        if ch == EOF: break
        if ch == <char>'\n':
            num_lines += 1

    print num_lines
    print num_factors
    return
    ptr_fr = fopen(fp.name, "rb")
    i = 0
    cdef int[::1] ptr_labels = numpy.empty((num_lines*num_factors,), dtype='int32')
    for i in xrange(num_lines):
        fscanf(ptr_fr,"%*s\t%i\t%i", &start, &stop)
        for j in xrange(num_factors):
            fscanf(ptr_fr,"\t%i", &ptr_labels[i*num_factors+j])
        i += 1
    return numpy.asarray(ptr_labels).reshape((num_lines, num_factors))
