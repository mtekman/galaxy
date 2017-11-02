#!/usr/bin/env python2
from bisect import bisect

class Converter:

    def __init__(self):
        self.out_lod     = "linkage.allegro_lod"
        self.out_haplo   = "linkage.allegro_haplo"
        self.out_descent = "linkage.allegro_descent"

        self.lod_map     = {}
        self.haplo_map   = {}
        self.descent_map = {}

        self.markers = {}
        self.marker_ordering = []


    @staticmethod
    def tokenizer(line):
        return line.splitlines()[0].split()
    
    # Override
    def extractLOD(self, file):
        return ""

    # Override
    def extractHaplo(self, file):
        return ""

    # Override
    def extractDescent(self, file):
        return ""
    
    def writeLOD(self):
        for variant in self.lod_map:
            # work

    def writeDescent(self):
        self.writeHaplo(True):

    def writeHaplo(self, descent=False):
        # write headers
        
        for indiv in self.haplo_map:
            # work

        



class Genehunter(Converter):


    def __init__(self, lod = None, haplo = None, descent = None):
        super(self)
        
        if lod is not None:
            self.extractLOD(lod)
            self.writeLOD()

        if haplo is not None:
            self.extractHaplo(haplo)
            self.writeHaplo()

        if descent is not None:
            self.extractDescent(descent)
            self.writeDescent()
            

    def extractLOD(self, file):

        pos_lod = {}  # pos → lod
        gt_line = ""
        
        with open(file,'r') as fio:

            # Jump to marker positions
            line = ""
            while not line.startswith("Current map ("):
                line = fio.readline()
                # iterate

            # Extract gts
            while not line.startswith("npl"):
                line = fio.readline()
                gt_line += line.splitlines()[0] + " "

            # Jump to LOD scoring
            while not line.startswith("position  LOD score    NPL score  p-value"):
                line = fio.readline()

            # Extract lod
            while not line.startswith(" "):
                gpos, lod, npl, pval, info  = map(float,Converter.tokenizer(line))

                if gpos not in pos_lod:
                    pos_lod[gpos] = lod
                else:
                    # update if new lod is larger
                    old_lod = pos_lod[gpos]
                    if lod > old_lod:
                        pos_lod[gpos] = lod


            fio.close()

        # Split marker pos pairs into map
        pos_pairs = gt_line.split()
        if len(pos_pairs) % 2 == 1:
            print >> sys.stderr, "even number of elements"
            exit(-1)
        
        pos_rsid = smap()
        cum_pos = 0
        ein = 0
        while ein < len(pos_pairs) - 1:
            marker = pos_pairs[ein]
            positn = float(pos_pairs[ein+1])

            pos_rsid[cum_pos] = marker

            cum_pos += positn
            ein += 2

        pos_rsid[cum_pos] = pos_pairs[-1]

        # Produce sorted list of gpos
        pos_rsid_array = sorted(pos_rsid.keys())
        pos_lod_array  = sorted(pos_lod.keys())

        keys_to_annotate = {}

        for marker_pos in pos_rsid_array:
            marker = pos_rsid[marker_pos]
            lm_ind = bisect_left(pos_lod_array, marker_pos)
            lm_pos = pos_lod_array[lm_ind - 1]

            if lm_pos not in keys_to_annotate:
                keys_to_annotate[lm_pos] = marker
            else:
                if marker_pos - 

                
            

        # Produce fam, loc, lod, marker
        ke
        

            
            
            
                
                
                
            
            
        


        
