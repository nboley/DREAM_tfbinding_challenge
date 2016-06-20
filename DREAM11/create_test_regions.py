import os, sys
import tempfile

from collections import namedtuple

import numpy as np

import h5py
import pybedtools

NarrowPeakData = namedtuple(
    'NarrowPeak', ['contig', 'start', 'stop', 'summit', 
                   'score', 'signalValue', 'pValue', 'qValue', 'idrValue', 'seq'])
NarrowPeakData.__new__.__defaults__ = (None,) * len(NarrowPeakData._fields)

class NarrowPeak(NarrowPeakData):
    @property
    def identifier(self):
        return "{0.contig}:{0.start}-{0.stop}_{0.summit}".format(self)
    @property    
    def pk_width(self):
        return self.stop - self.start

class Regions(object):
    @property
    def regions_bed(self):
        if not hasattr(self, '_regions_bed'):
            self._regions_bed = save_regions_array_to_bedtool(self.regions)
        return self._regions_bed

    @property
    def regions_with_flank_bed(self):
        if not hasattr(self, '_regions_with_flank_bed'):
            self._regions_with_flank_bed = save_regions_array_to_bedtool(
                self.regions_with_flank)
        return self._regions_with_flank_bed

    def save(self, ofname):        
        with h5py.File(ofname, "w") as f:
            f['regions'] = self.regions
            f['regions_with_flank'] = self.regions_with_flank
        return

    @classmethod
    def load(cls, fname):
        f = h5py.File(fname, 'r')
        return cls(f['regions'], f['regions_with_flank'])

    def __init__(self, regions, regions_with_flank=None):
        assert isinstance(regions, (np.ndarray, h5py._hl.dataset.Dataset))
        self.regions = regions
        if regions_with_flank is None:
            self.regions_with_flank = self.regions
        else:
            self.regions_with_flank = regions_with_flank
        return

def save_regions_array_to_bedtool(regions):
    bed_fp = tempfile.NamedTemporaryFile("w", dir="./", delete=False, suffix=".bed")
    np.savetxt(bed_fp, regions, fmt="%s\t%i\t%i")
    bed_fp.flush()
    return pybedtools.BedTool(bed_fp.name)
    
def iter_expanded_regions(regions, min_size):
    """Expand regions to cover at least min_size bases.

    """
    for region in regions:
        start, stop = region[1], region[2]
        pk_width = stop-start
        if pk_width < min_size:
            start = max(0, start-(pk_width-min_size)/2)
            stop = start+min_size
        yield NarrowPeak(region[0], start, stop)
    return

def iter_rounded_regions(regions, min_block_size):
    """Expand regions so that region_size.(start/stop)%block_size == 0.

    This is used to ensure that tiled regions dont intersect. 
    """
    for region in regions:
        start, stop = region[1], region[2]
        start = min_block_size*(start//min_block_size)
        if stop%min_block_size > 0:
            stop += (min_block_size - stop%min_block_size)
        yield NarrowPeak(region[0], start, stop)
    return

def build_genomic_regions_array(regions):
    interval_types = ('S64', 'i4', 'i4')
    return np.array(
        [(str(x[0]), int(x[1]), int(x[2])) for x in regions], 
        dtype=zip(('contig', 'start', 'stop'), interval_types)
    )

def build_genomic_regions_array_from_bed(regions):
    return build_genomic_regions_array(regions)

def build_tiled_regions(regions, core_size, flank_size, offset_size):
    """Expand and iterate tiled regions.
    
    """
    assert core_size%offset_size == 0
    new_region_size = core_size + 2*flank_size
    assert new_region_size%offset_size == 0
    regions_str = "\n".join(
        "\t".join(map(str, x[:3])) 
        for x in iter_rounded_regions(
            iter_expanded_regions(regions, new_region_size),
            offset_size
        )
    )
    regions_bed = pybedtools.BedTool(regions_str, from_string=True)
    regions_bed = regions_bed.sort()
    regions_bed = regions_bed.merge()
    tiled_regions = []
    tiled_regions_with_flank = []
    for region in regions_bed:
        region_size = region.stop-region.start
        assert region_size%offset_size == 0
        num_offsets = region_size//offset_size
        for i in xrange(num_offsets):
            start = region.start-flank_size+i*offset_size
            if start < 0: continue
            stop = region.start+core_size+flank_size+i*offset_size
            tiled_regions_with_flank.append(NarrowPeak(region[0], start, stop))
            tiled_regions.append(
                NarrowPeak(region[0], start+flank_size, stop-flank_size))

    regions_str = "\n".join(
        ["\t".join(map(str, x[:3])) for x in tiled_regions])
    regions_bed = pybedtools.BedTool(regions_str, from_string=True)
    regions_bed = regions_bed.sort()

    regions_w_flank_str = "\n".join(
        ["\t".join(map(str, x[:3])) for x in tiled_regions_with_flank])
    regions_w_flank_bed = pybedtools.BedTool(regions_w_flank_str, from_string=True)
    regions_w_flank_bed = regions_w_flank_bed.sort()
    return Regions(
        build_genomic_regions_array_from_bed(regions_bed), 
        build_genomic_regions_array_from_bed(regions_w_flank_bed)
    )

def create_genomic_regions(core_size=200, flank_size=400, offset_size=50):
    # primary contigs and their lengths for hg19 
    # dropped for consistency chrY	59373566
    contig_regions = tuple([
        (x.split()[0],
         int(x.split()[1])+flank_size+core_size,
         int(x.split()[2])-flank_size-core_size)
        for x in 
    """
chr1	0	249250621
chr2	0	243199373
chr3	0	198022430
chr4	0	191154276
chr5	0	180915260
chr6	0	171115067
chr7	0	159138663
chr8	0	146364022
chr9	0	141213431
chr10	0	135534747
chr11	0	135006516
chr12	0	133851895
chr13	0	115169878
chr14	0	107349540
chr15	0	102531392
chr16	0	90354753
chr17	0	81195210
chr18	0	78077248
chr20	0	63025520
chr19	0	59128983
chr21	0	48129895
chr22	0	51304566
chrX	0	155270560
    """.strip().split("\n")
    ])
    #contig_regions = contig_regions[-2:]
    cached_fname = "all_regions.%s.h5" % abs(
        hash((contig_regions, core_size, flank_size, offset_size)))
    try:
        return Regions.load(cached_fname)
    except IOError:
        regions = build_tiled_regions(
            contig_regions, core_size, flank_size, offset_size)
        regions.save(cached_fname)
        return regions

def main():
    regions = create_genomic_regions()
    print regions
    with open("all_regions.bed", "w") as ofp:
        np.savetxt(ofp, regions.regions, fmt="%s\t%i\t%i")

if __name__ == '__main__':
    main()
