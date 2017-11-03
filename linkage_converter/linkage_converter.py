#!/usr/bin/env python2
from bisect import bisect_left
import sys

class Converter:
    """
    Superclass for converting formats
    """
    haplo_individual_buffer = "%5d %5d %5d %5d %2d %2d"
    lod_line_buffer = "%8s %8s %8s %8s %-s"

    def __init__(self, pedin, mapin):
        self.out_lod = "linkage.allegro_lod"
        self.out_haplo = "linkage.allegro_haplo"
        self.out_descent = "linkage.allegro_descent"

        self.haplo_map = {}  # fam_id → indiv_id → (all1,all2)
        self.pos_marker = {} # genpos → marker

        self.__populatePedigree(pedin)
        self.__populateMarkerMap(mapin)


    def __populatePedigree(input_ped):
        with open(input_ped, "r") as pio:
            for line in pio:
                tokens = map(int, line.split())
                f_id, p_id, father, mother, gender, affect = tokens[:6]

                if f_id not in self.pedigree:
                    self.pedigree[f_id] = {}
                if p_id not in self.pedigree:
                    self.pedigree[f_id][p_id] = (father, mother, gender, affect)
                else:
                    print >> sys.stderr, "Duplicate individual:", f_id, p_id
                    exit(-1)


    def __populateMarkerMap(mapin):
        with open(mapin, "r") as mio:
            mio.readline() # chomp header

            for line in mio:
                chrom, gpos, marker, ppos, nr = Converter.tokenizer(line)
                self.pos_marker[float(gpos)] = marker.strip()


    @staticmethod
    def tokenizer(line):
        return line.splitlines()[0].split()

    # Override
    def extractLOD(self, file1):
        return ""

    # Override
    def extractHaplo(self, file1, extrafile=None):
        return ""

    # Override
    def extractDescent(self, file1):
        return ""

    def writeLOD(self):

        out_lod = open(self.out_lod, "w")

        header = Converter.lod_line_buffer % (
            "location", "LOD", "alpha", "HLOD", "marker"
        )
        print(header, file=out_lod)

        for pos, lod, alpha, hlod, marker in self.lod_array:
            out_line = Converter.lod_line_buffer % (
                "%.4f" % pos,
                "%.4f" % lod,
                "%.4f" % alpha,
                "%.4f" % hlod,
                marker
            )
            print(out_line, file=out_lod)
        out_lod.close()



    def writeDescent(self):
        self.writeHaplo(True)

    def writeHaplo(self, descent=False):

        out_haplo = open(self.out_haplo, "w")

        headers = self.__generateHeaders()
        print(headers, file=out_haplo)

        for fam_id in self.haplo_map:
            for indiv_id in self.haplo_map[fam_id]:
                alleles = self.haplo_map[fam_id][indiv_id]
                ped_data = self.pedigree[fam_id][indiv_id]
                #father, mother, gender, affect = ped_data

                indiv_data = Converter.haplo_individual_buffer % (
                    fam_id, indiv_id, *ped_data
                )

                print(indiv_data, alleles[0], file=out_haplo)
                print(indiv_data, alleles[0], file=out_haplo)
        out_haplo.close()



class Simwalk(Converter):

    def __init__(self, pedigree, mapfile, hef):
        super(self, pedigree, mapfile)


    def extractLOD(

        



class Genehunter(Converter):

    def __init__(self, pedigree, mapfile, lod=None, haplo=None):
        super(self, pedigree, mapfile)

        if lod is not None:
            self.extractLOD(lod)
            self.writeLOD()

        if haplo is not None:
            self.extractHaplo(haplo)
            self.writeHaplo()


    #def __processMarkers(gt_line, pos_lod):
    def __processMarkers(self, pos_lod):
        # # Split marker pos pairs into map
        # pos_pairs = gt_line.split()
        # if len(pos_pairs) % 2 == 1:
        #     print >> sys.stderr, "even number of elements"
        #     exit(-1)

        # pos_rsid = {}
        # cum_pos = 0
        # ein = 0
        # while ein < len(pos_pairs) - 1:
        #     marker = pos_pairs[ein]
        #     positn = float(pos_pairs[ein+1])

        #     pos_rsid[cum_pos] = marker

        #     cum_pos += positn
        #     ein += 2

        # pos_rsid[cum_pos] = pos_pairs[-1]

        # Produce sorted list of gpos
        #pos_rsid_array = sorted(pos_rsid.keys())
        pos_rsid_array = sorted(self.pos_marker.keys())                                
        pos_lod_array = sorted(pos_lod.keys())

        keys_to_annotate = {}

        # bisect and find closest marker
        for marker_pos in pos_rsid_array:
            marker = self.pos_marker[marker_pos]
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


    def extractHaplo(self, hapfile):
        with open(hapfile, 'r') as hio:
            line = "   "
            fam_id = -1

            try:
                while True:
                    while not line.startswith("*****"):
                        line = hio.readline()
                        fam_id = int(line.split()[1])

                        if fam_id not in self.indiv_map:
                            self.haplo_map[fam_id] = {}

                    while not line.startswith(" "):
                        line = hio.readline()
                        indiv_all1 = map(int, line.split())
                        line = hio.readline()
                        indiv_all2 = map(int, line.split())

                        indiv_id = indiv_all1[0]

                        if len(indiv_all1) != len(indiv_all2) + 4:
                            print >> sys.stderr, "Allele mismatch for indiv", fam_id, indiv_id
                            exit(-1)

                        self.haplo_map[fam_id][indiv_id] = (indiv_all1[4:], indiv_all2)
            except IOError:
                hio.close()


    def extractLOD(self, file):

        pos_lod = {}  # pos → lod
        gt_line = ""
        num_peds = 0

        with open(file, 'r') as fio:

            # Jump to marker positions
            line = ""
            while not line.startswith("Current map ("):
                line = fio.readline()
                # iterate

            # # Extract gts
            # while not line.startswith("npl"):
            #     line = fio.readline()
            #     gt_line += line.splitlines()[0] + " "

            # Get num pedigrees
            while not line.startswith("Totalling pedigrees"):
                line = fio.readline()
                # iterate

            # Jump to LOD scoring
            while not line.startswith("position LOD_score  (alpha, HLOD)"):
                line = fio.readline()

            # Extract lod
            while not line.startswith(" "):
                gpos, lod, alpha, hlod, npl, pval, info = Converter.tokenizer(line)
                gpos = float(gpos)

                if lod == "-INFINITY":
                    lod = -10000
                lod = float(lod)

                alpha = float(alpha.split("(")[-1])
                hlod = float(hlod.split(")")[0])

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


        #postns_with_closest_markers, sorted_poslod = self.__processMarkers(gt_line, pos_lod)
        postns_with_closest_markers, sorted_poslod = self.__processMarkers(pos_lod)

        # update map with marker info
        for pos in sorted_poslod:
            marker = ""
            if pos in postns_with_closest_markers:
                marker = postns_with_closest_markers[pos]

            lod, alpha, hlod, npl, pval = pos_lod[pos]

            self.lod_array.append((pos, lod, alpha, hlod, marker))
