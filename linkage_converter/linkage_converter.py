#!/usr/bin/env python2
from bisect import bisect

class Converter:

    def __init__(self):
        self.out_lod     = "linkage.allegro_lod"
        self.out_haplo   = "linkage.allegro_haplo"
        self.out_descent = "linkage.allegro_descent"

        self.lod_array     = {}
        self.haplo_map   = {}
        self.descent_map = {}

        self.markers = {}
        self.marker_ordering = []


    @staticmethod
    def tokenizer(line):
        return line.splitlines()[0].split()
    
    # Override
    def extractLOD(self, file1):
        return ""

    # Override
    def extractHaplo(self, file1, extrafile = None):
        return ""

    # Override
    def extractDescent(self, file1):
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


    def __init__(self, lod = None, haplo = None, descent = None, mapfile = None):
        super(self)
        
        if lod is not None:
            self.extractLOD(lod)
            self.writeLOD()

        if haplo is not None:
            if mapfile is None:
                print >> sys.stderr, "map file required"
                exit(-1)
                
            self.extractHaplo(haplo, mapfile)
            self.writeHaplo()

        if descent is not None:
            self.extractDescent(descent)
            self.writeDescent()


    def __processMarkers(gt_line, pos_lod):
        # Split marker pos pairs into map
        pos_pairs = gt_line.split()
        if len(pos_pairs) % 2 == 1:
            print >> sys.stderr, "even number of elements"
            exit(-1)
        
        pos_rsid = {}
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

        # bisect and find closest marker
        for marker_pos in pos_rsid_array:
            marker = pos_rsid[marker_pos]
            lm_ind = bisect_left(pos_lod_array, marker_pos)
            lm_pos = pos_lod_array[lm_ind - 1]

            diff1 = marker_pos - lm_pos
            diff1 = diff1 if diff1 > 0 else -diff1

            if diff1 > 0.1:
                continue

            if lm_pos not in keys_to_annotate:
                keys_to_annotate[lm_pos] = marker
            else:
                diff2 = marker_pos - keys_to_annotate[lm_pos]             
                diff2 = diff2 if diff2 > 0 else -diff2

                if diff1 < diff2:
                    # update with closer marker
                    keys_to_annotate[lm_pos] = marker


        return keys_to_annotate, pos_lod_array


    def extractHaplo(self, hapfile, mapfile):
        with open(haplfile, 'r') as hio:
            hio.readline()

            # read data in pairs
            line = "   "
            while not line.startswith(" "):
                
            
            


    def extractLOD(self, file):

        pos_lod = {}  # pos → lod
        gt_line = ""
        num_peds = 0
        
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

            # Get num pedigrees
            while not line.startswith("Totalling pedigrees"):
                line = fio.readline()
                # iterate

            # Jump to LOD scoring
            while not line.startswith("position LOD_score  (alpha, HLOD)"):
                line = fio.readline()

            # Extract lod
            while not line.startswith(" "):
                gpos, lod, alpha, hlod, npl, pval, info  = Converter.tokenizer(line))
                gpos = float(gpos)

                if lod == "-INFINITY":
                    lod = -10000
                lod = float(lod)

                alpha = float(alpha.split("(")[-1])
                hlod  = float(hlod.split(")")[0])

                npl = float(npl)
                pval = float(pval)

                # insert
                if gpos not in pos_lod:
                    pos_lod[gpos] = [lod, alpha, hlod, npl, pval]
                else:
                    # update if new lod is larger
                    old_lod = pos_lod[gpos][0]
                    if lod > old_lod:
                        pos_lod[gpos][0] = lod

            fio.close()


        postns_with_closest_markers, sorted_poslod = self.__processMarkers(gt_line, pos_lod)

        # update map with marker info
        for pos in sorted_poslod:
            marker = ""
            if pos in postns_with_closest_markers:
                marker = postns_with_closest_markers[pos]

            lod, alpha, hlod, npl, pval = pos_lod[pos]
            
            self.lod_array.append( (pos, lod, alpha, hlod, marker) )
            
