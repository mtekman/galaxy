"""
rgenetics datatypes
Use at your peril
Ross Lazarus
for the rgenetics and galaxy projects

genome graphs datatypes derived from Interval datatypes
genome graphs datasets have a header row with appropriate columnames
The first column is always the marker - eg columname = rs, first row= rs12345 if the rows are snps
subsequent row values are all numeric ! Will fail if any non numeric (eg '+' or 'NA') values
ross lazarus for rgenetics
august 20 2007
"""
import logging
import os
import re
import sys
from cgi import escape

from six.moves.urllib.parse import quote_plus

from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.text import Html
from galaxy.datatypes.data import Text
from galaxy.util import nice_size
from galaxy.web import url_for

gal_Log = logging.getLogger(__name__)
verbose = False

# https://genome.ucsc.edu/goldenpath/help/hgGenomeHelp.html
VALID_GENOME_GRAPH_MARKERS = re.compile('^(chr.*|RH.*|rs.*|SNP_.*|CN.*|A_.*)')
VALID_GENOTYPES_LINE = re.compile('^([a-zA-Z0-9]+)(\\s([0-9]{2}|[A-Z]{2}|NC|\?\?))+\\s*$')


class GenomeGraphs(Tabular):
    """
    Tab delimited data containing a marker id and any number of numeric values
    """

    MetadataElement(name="markerCol", default=1, desc="Marker ID column", param=metadata.ColumnParameter)
    MetadataElement(name="columns", default=3, desc="Number of columns", readonly=True)
    MetadataElement(name="column_types", default=[], desc="Column types", readonly=True, visible=False)
    file_ext = 'gg'

    def __init__(self, **kwd):
        """
        Initialize gg datatype, by adding UCSC display apps
        """
        Tabular.__init__(self, **kwd)
        self.add_display_app('ucsc', 'Genome Graph', 'as_ucsc_display_file', 'ucsc_links')

    def set_meta(self, dataset, **kwd):
        Tabular.set_meta(self, dataset, **kwd)
        dataset.metadata.markerCol = 1
        header = open(dataset.file_name, 'r').readlines()[0].strip().split('\t')
        dataset.metadata.columns = len(header)
        t = ['numeric' for x in header]
        t[0] = 'string'
        dataset.metadata.column_types = t
        return True

    def as_ucsc_display_file(self, dataset, **kwd):
        """
        Returns file
        """
        return open(dataset.file_name, 'r')

    def ucsc_links(self, dataset, type, app, base_url):
        """
        from the ever-helpful angie hinrichs angie@soe.ucsc.edu
        a genome graphs call looks like this

        http://genome.ucsc.edu/cgi-bin/hgGenome?clade=mammal&org=Human&db=hg18&hgGenome_dataSetName=dname
        &hgGenome_dataSetDescription=test&hgGenome_formatType=best%20guess&hgGenome_markerType=best%20guess
        &hgGenome_columnLabels=best%20guess&hgGenome_maxVal=&hgGenome_labelVals=
        &hgGenome_maxGapToFill=25000000&hgGenome_uploadFile=http://galaxy.esphealth.org/datasets/333/display/index
        &hgGenome_doSubmitUpload=submit

        Galaxy gives this for an interval file

        http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr1:1-1000&hgt.customText=
        http%3A%2F%2Fgalaxy.esphealth.org%2Fdisplay_as%3Fid%3D339%26display_app%3Ducsc

        """
        ret_val = []
        if not dataset.dbkey:
            dataset.dbkey = 'hg18'  # punt!
        if dataset.has_data():
            for site_name, site_url in app.datatypes_registry.get_legacy_sites_by_build('ucsc', dataset.dbkey):
                if site_name in app.datatypes_registry.get_display_sites('ucsc'):
                    site_url = site_url.replace('/hgTracks?', '/hgGenome?')  # for genome graphs
                    internal_url = "%s" % url_for(controller='dataset',
                                                  dataset_id=dataset.id,
                                                  action='display_at',
                                                  filename='ucsc_' + site_name)
                    display_url = "%s%s/display_as?id=%i&display_app=%s&authz_method=display_at" % (base_url, url_for(controller='root'), dataset.id, type)
                    display_url = quote_plus(display_url)
                    # was display_url = quote_plus( "%s/display_as?id=%i&display_app=%s" % (base_url, dataset.id, type) )
                    # redirect_url = quote_plus( "%sdb=%s&position=%s:%s-%s&hgt.customText=%%s" % (site_url, dataset.dbkey, chrom, start, stop) )
                    sl = ["%sdb=%s" % (site_url, dataset.dbkey), ]
                    # sl.append("&hgt.customText=%s")
                    sl.append("&hgGenome_dataSetName=%s&hgGenome_dataSetDescription=%s" % (dataset.name, 'GalaxyGG_data'))
                    sl.append("&hgGenome_formatType=best guess&hgGenome_markerType=best guess")
                    sl.append("&hgGenome_columnLabels=first row&hgGenome_maxVal=&hgGenome_labelVals=")
                    sl.append("&hgGenome_doSubmitUpload=submit")
                    sl.append("&hgGenome_maxGapToFill=25000000&hgGenome_uploadFile=%s" % display_url)
                    s = ''.join(sl)
                    s = quote_plus(s)
                    redirect_url = s
                    link = '%s?redirect_url=%s&display_url=%s' % (internal_url, redirect_url, display_url)
                    ret_val.append((site_name, link))
        return ret_val

    def make_html_table(self, dataset, skipchars=[]):
        """
        Create HTML table, used for displaying peek
        """
        out = ['<table cellspacing="0" cellpadding="3">']
        try:
            with open(dataset.file_name, 'r') as f:
                d = f.readlines()[:5]
            if len(d) == 0:
                out = "Cannot find anything to parse in %s" % dataset.name
                return out
            hasheader = 0
            try:
                ['%f' % x for x in d[0][1:]]  # first is name - see if starts all numerics
            except Exception:
                hasheader = 1
            # Generate column header
            out.append('<tr>')
            if hasheader:
                for i, name in enumerate(d[0].split()):
                    out.append('<th>%s.%s</th>' % (i + 1, name))
                d.pop(0)
                out.append('</tr>')
            for row in d:
                out.append('<tr>')
                out.append(''.join(['<td>%s</td>' % x for x in row.split()]))
                out.append('</tr>')
            out.append('</table>')
            out = "".join(out)
        except Exception as exc:
            out = "Can't create peek %s" % exc
        return out

    def validate(self, dataset):
        """
        Validate a gg file - all numeric after header row
        """
        errors = list()
        with open(dataset.file_name, "r") as infile:
            next(infile)  # header
            for i, row in enumerate(infile):
                ll = row.strip().split('\t')[1:]  # first is alpha feature identifier
                badvals = []
                for j, x in enumerate(ll):
                    try:
                        x = float(x)
                    except Exception:
                        badvals.append('col%d:%s' % (j + 1, x))
        if len(badvals) > 0:
            errors.append('row %d, %s' % (' '.join(badvals)))
            return errors

    def sniff(self, filename):
        """
        Determines whether the file is in gg format

        >>> from galaxy.datatypes.sniff import get_test_fname
        >>> fname = get_test_fname( 'test_space.txt' )
        >>> GenomeGraphs().sniff( fname )
        False
        >>> fname = get_test_fname( '1.gg' )
        >>> GenomeGraphs().sniff( fname )
        True
        """
        with open(filename, 'r') as f:
            buf = f.read(1024)

        rows = [l.split() for l in buf.splitlines()[1:4]]  # break on lines and drop header, small sample

        if len(rows) < 1:
            return False

        for row in rows:
            if len(row) < 2:
                # Must actually have a marker and at least one numeric value
                return False
            first_val = row[0]
            if not VALID_GENOME_GRAPH_MARKERS.match(first_val):
                return False
            rest_row = row[1:]
            try:
                [float(x) for x in rest_row]  # first col has been removed
            except ValueError:
                return False
        return True

    def get_mime(self):
        """Returns the mime type of the datatype"""
        return 'application/vnd.ms-excel'


class rgTabList(Tabular):
    """
    for sampleid and for featureid lists of exclusions or inclusions in the clean tool
    featureid subsets on statistical criteria -> specialized display such as gg
    """
    file_ext = "rgTList"

    def __init__(self, **kwd):
        """
        Initialize featurelistt datatype
        """
        Tabular.__init__(self, **kwd)
        self.column_names = []

    def display_peek(self, dataset):
        """Returns formated html of peek"""
        return self.make_html_table(dataset, column_names=self.column_names)

    def get_mime(self):
        """Returns the mime type of the datatype"""
        return 'text/html'


class rgSampleList(rgTabList):
    """
    for sampleid exclusions or inclusions in the clean tool
    output from QC eg excess het, gender error, ibd pair member,eigen outlier,excess mendel errors,...
    since they can be uploaded, should be flexible
    but they are persistent at least
    same infrastructure for expression?
    """
    file_ext = "rgSList"

    def __init__(self, **kwd):
        """
        Initialize samplelist datatype
        """
        rgTabList.__init__(self, **kwd)
        self.column_names[0] = 'FID'
        self.column_names[1] = 'IID'
        # this is what Plink wants as at 2009

    def sniff(self, filename):
        with open(filename, "r") as infile:
            header = next(infile)  # header
        if header[0] == 'FID' and header[1] == 'IID':
            return True
        else:
            return False


class rgFeatureList(rgTabList):
    """
    for featureid lists of exclusions or inclusions in the clean tool
    output from QC eg low maf, high missingness, bad hwe in controls, excess mendel errors,...
    featureid subsets on statistical criteria -> specialized display such as gg
    same infrastructure for expression?
    """
    file_ext = "rgFList"

    def __init__(self, **kwd):
        """Initialize featurelist datatype"""
        rgTabList.__init__(self, **kwd)
        for i, s in enumerate(['#FeatureId', 'Chr', 'Genpos', 'Mappos']):
            self.column_names[i] = s


class Rgenetics(Html):
    """
    base class to use for rgenetics datatypes
    derived from html - composite datatype elements
    stored in extra files path
    """

    MetadataElement(name="base_name", desc="base name for all transformed versions of this genetic dataset", default='RgeneticsData',
                    readonly=True, set_in_upload=True)

    composite_type = 'auto_primary_file'
    allow_datatype_change = False
    file_ext = 'rgenetics'

    def generate_primary_file(self, dataset=None):
        rval = ['<html><head><title>Rgenetics Galaxy Composite Dataset </title></head><p/>']
        rval.append('<div>This composite dataset is composed of the following files:<p/><ul>')
        for composite_name, composite_file in self.get_composite_files(dataset=dataset).items():
            fn = composite_name
            opt_text = ''
            if composite_file.optional:
                opt_text = ' (optional)'
            if composite_file.get('description'):
                rval.append('<li><a href="%s" type="application/binary">%s (%s)</a>%s</li>' % (fn, fn, composite_file.get('description'), opt_text))
            else:
                rval.append('<li><a href="%s" type="application/binary">%s</a>%s</li>' % (fn, fn, opt_text))
        rval.append('</ul></div></html>')
        return "\n".join(rval)

    def regenerate_primary_file(self, dataset):
        """
        cannot do this until we are setting metadata
        """
        efp = dataset.extra_files_path
        flist = os.listdir(efp)
        rval = ['<html><head><title>Files for Composite Dataset %s</title></head><body><p/>Composite %s contains:<p/><ul>' % (dataset.name, dataset.name)]
        for i, fname in enumerate(flist):
            sfname = os.path.split(fname)[-1]
            f, e = os.path.splitext(fname)
            rval.append('<li><a href="%s">%s</a></li>' % (sfname, sfname))
        rval.append('</ul></body></html>')
        with open(dataset.file_name, 'w') as f:
            f.write("\n".join(rval))
            f.write('\n')

    def get_mime(self):
        """Returns the mime type of the datatype"""
        return 'text/html'

    def set_meta(self, dataset, **kwd):
        """
        for lped/pbed eg

        """
        Html.set_meta(self, dataset, **kwd)
        if not kwd.get('overwrite'):
            if verbose:
                gal_Log.debug('@@@ rgenetics set_meta called with overwrite = False')
            return True
        try:
            efp = dataset.extra_files_path
        except Exception:
            if verbose:
                gal_Log.debug('@@@rgenetics set_meta failed %s - dataset %s has no efp ?' % (sys.exc_info()[0], dataset.name))
            return False
        try:
            flist = os.listdir(efp)
        except Exception:
            if verbose:
                gal_Log.debug('@@@rgenetics set_meta failed %s - dataset %s has no efp ?' % (sys.exc_info()[0], dataset.name))
            return False
        if len(flist) == 0:
            if verbose:
                gal_Log.debug('@@@rgenetics set_meta failed - %s efp %s is empty?' % (dataset.name, efp))
            return False
        self.regenerate_primary_file(dataset)
        if not dataset.info:
            dataset.info = 'Galaxy genotype datatype object'
        if not dataset.blurb:
            dataset.blurb = 'Composite file - Rgenetics Galaxy toolkit'
        return True


class SNPMatrix(Rgenetics):
    """
    BioC SNPMatrix Rgenetics data collections
    """
    file_ext = "snpmatrix"

    def set_peek(self, dataset, **kwd):
        if not dataset.dataset.purged:
            dataset.peek = "Binary RGenetics file"
            dataset.blurb = nice_size(dataset.get_size())
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def sniff(self, filename):
        """ need to check the file header hex code
        """
        infile = open(filename, "b")
        head = infile.read(16)
        head = [hex(x) for x in head]
        if head != '':
            return False
        else:
            return True


class Lped(Rgenetics):
    """
    linkage pedigree (ped,map) Rgenetics data collections
    """
    file_ext = "lped"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.ped',
                                description='Pedigree File',
                                substitute_name_with_metadata='base_name',
                                is_binary=False)
        self.add_composite_file('%s.map',
                                description='Map File',
                                substitute_name_with_metadata='base_name',
                                is_binary=False)


class Pphe(Rgenetics):
    """
    Plink phenotype file - header must have FID\tIID... Rgenetics data collections
    """
    file_ext = "pphe"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.pphe',
                                description='Plink Phenotype File',
                                substitute_name_with_metadata='base_name',
                                is_binary=False)


class Fphe(Rgenetics):
    """
    fbat pedigree file - mad format with ! as first char on header row
    Rgenetics data collections
    """
    file_ext = "fphe"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.fphe',
                                description='FBAT Phenotype File',
                                substitute_name_with_metadata='base_name')


class Phe(Rgenetics):
    """
    Phenotype file
    """
    file_ext = "phe"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.phe',
                                description='Phenotype File',
                                substitute_name_with_metadata='base_name',
                                is_binary=False)


class Fped(Rgenetics):
    """
    FBAT pedigree format - single file, map is header row of rs numbers. Strange.
    Rgenetics data collections
    """
    file_ext = "fped"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.fped', description='FBAT format pedfile',
                                substitute_name_with_metadata='base_name',
                                is_binary=False)


class Pbed(Rgenetics):
    """
    Plink Binary compressed 2bit/geno Rgenetics data collections
    """
    file_ext = "pbed"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.bim', substitute_name_with_metadata='base_name', is_binary=False)
        self.add_composite_file('%s.bed', substitute_name_with_metadata='base_name', is_binary=True)
        self.add_composite_file('%s.fam', substitute_name_with_metadata='base_name', is_binary=False)


class ldIndep(Rgenetics):
    """
    LD (a good measure of redundancy of information) depleted Plink Binary compressed 2bit/geno
    This is really a plink binary, but some tools work better with less redundancy so are constrained to
    these files
    """
    file_ext = "ldreduced"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.bim', substitute_name_with_metadata='base_name', is_binary=False)
        self.add_composite_file('%s.bed', substitute_name_with_metadata='base_name', is_binary=True)
        self.add_composite_file('%s.fam', substitute_name_with_metadata='base_name', is_binary=False)


class Eigenstratgeno(Rgenetics):
    """
    Eigenstrat format - may be able to get rid of this
    if we move to shellfish
    Rgenetics data collections
    """
    file_ext = "eigenstratgeno"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.eigenstratgeno', substitute_name_with_metadata='base_name', is_binary=False)
        self.add_composite_file('%s.ind', substitute_name_with_metadata='base_name', is_binary=False)
        self.add_composite_file('%s.map', substitute_name_with_metadata='base_name', is_binary=False)


class Eigenstratpca(Rgenetics):
    """
    Eigenstrat PCA file for case control adjustment
    Rgenetics data collections
    """
    file_ext = "eigenstratpca"

    def __init__(self, **kwd):
        Rgenetics.__init__(self, **kwd)
        self.add_composite_file('%s.eigenstratpca',
                                description='Eigenstrat PCA file', substitute_name_with_metadata='base_name')


class Snptest(Rgenetics):
    """
    BioC snptest Rgenetics data collections
    """
    file_ext = "snptest"


class Pheno(Tabular):
    """
    base class for pheno files
    """
    file_ext = 'pheno'


class RexpBase(Html):
    """
    base class for BioC data structures in Galaxy
    must be constructed with the pheno data in place since that
    goes into the metadata for each instance
    """
    MetadataElement(name="columns", default=0, desc="Number of columns", visible=True)
    MetadataElement(name="column_names", default=[], desc="Column names", visible=True)
    MetadataElement(name="pheCols", default=[], desc="Select list for potentially interesting variables", visible=True)
    MetadataElement(name="base_name",
                    desc="base name for all transformed versions of this expression dataset", default='rexpression', set_in_upload=True)
    MetadataElement(name="pheno_path", desc="Path to phenotype data for this experiment", default="rexpression.pheno", visible=True)
    file_ext = 'rexpbase'
    html_table = None
    is_binary = True
    composite_type = 'auto_primary_file'
    allow_datatype_change = False

    def __init__(self, **kwd):
        Html.__init__(self, **kwd)
        self.add_composite_file('%s.pheno', description='Phenodata tab text file',
                                substitute_name_with_metadata='base_name', is_binary=False)

    def generate_primary_file(self, dataset=None):
        """
        This is called only at upload to write the html file
        cannot rename the datasets here - they come with the default unfortunately
        """
        return '<html><head></head><body>AutoGenerated Primary File for Composite Dataset</body></html>'

    def get_mime(self):
        """Returns the mime type of the datatype"""
        return 'text/html'

    def get_phecols(self, phenolist=[], maxConc=20):
        """
        sept 2009: cannot use whitespace to split - make a more complex structure here
        and adjust the methods that rely on this structure
        return interesting phenotype column names for an rexpression eset or affybatch
        to use in array subsetting and so on. Returns a data structure for a
        dynamic Galaxy select parameter.
        A column with only 1 value doesn't change, so is not interesting for
        analysis. A column with a different value in every row is equivalent to a unique
        identifier so is also not interesting for anova or limma analysis - both these
        are removed after the concordance (count of unique terms) is constructed for each
        column. Then a complication - each remaining pair of columns is tested for
        redundancy - if two columns are always paired, then only one is needed :)
        """
        for nrows, row in enumerate(phenolist):  # construct concordance
            if len(row.strip()) == 0:
                break
            row = row.strip().split('\t')
            if nrows == 0:  # set up from header
                head = row
                totcols = len(row)
                concordance = [{} for x in head]  # list of dicts
            else:
                for col, code in enumerate(row):  # keep column order correct
                    if col >= totcols:
                        gal_Log.warning('### get_phecols error in pheno file - row %d col %d (%s) longer than header %s' % (nrows, col, row, head))
                    else:
                        concordance[col].setdefault(code, 0)  # first one is zero
                        concordance[col][code] += 1
        useCols = []
        useConc = []  # columns of interest to keep
        nrows = len(phenolist)
        nrows -= 1  # drop head from count
        for c, conc in enumerate(concordance):  # c is column number
            if (len(conc) > 1) and (len(conc) < min(nrows, maxConc)):  # not all same and not all different!!
                useConc.append(conc)  # keep concordance
                useCols.append(c)  # keep column
        nuse = len(useCols)
        # now to check for pairs of concordant columns - drop one of these.
        delme = []
        p = phenolist[1:]  # drop header
        plist = [x.strip().split('\t') for x in p]  # list of lists
        phe = [[x[i] for i in useCols] for x in plist if len(x) >= totcols]  # strip unused data
        for i in range(0, (nuse - 1)):  # for each interesting column
            for j in range(i + 1, nuse):
                kdict = {}
                for row in phe:  # row is a list of lists
                    k = '%s%s' % (row[i], row[j])  # composite key
                    kdict[k] = k
                if (len(kdict.keys()) == len(concordance[useCols[j]])):  # i and j are always matched
                    delme.append(j)
        delme = list(set(delme))  # remove dupes
        listCol = []
        delme.sort()
        delme.reverse()  # must delete from far end!
        for i in delme:
            del useConc[i]  # get rid of concordance
            del useCols[i]  # and usecols entry
        for i, conc in enumerate(useConc):  # these are all unique columns for the design matrix
            ccounts = sorted((conc.get(code, 0), code) for code in conc.keys())  # decorate
            cc = [(x[1], x[0]) for x in ccounts]  # list of code count tuples
            codeDetails = (head[useCols[i]], cc)  # ('foo',[('a',3),('b',11),..])
            listCol.append(codeDetails)
        if len(listCol) > 0:
            res = listCol
            # metadata.pheCols becomes [('bar;22,zot;113','foo'), ...]
        else:
            res = [('no usable phenotype columns found', [('?', 0), ]), ]
        return res

    def get_pheno(self, dataset):
        """
        expects a .pheno file in the extra_files_dir - ugh
        note that R is wierd and adds the row.name in
        the header so the columns are all wrong - unless you tell it not to.
        A file can be written as
        write.table(file='foo.pheno',pData(foo),sep='\t',quote=F,row.names=F)
        """
        p = open(dataset.metadata.pheno_path, 'r').readlines()
        if len(p) > 0:  # should only need to fix an R pheno file once
            head = p[0].strip().split('\t')
            line1 = p[1].strip().split('\t')
            if len(head) < len(line1):
                head.insert(0, 'ChipFileName')  # fix R write.table b0rken-ness
                p[0] = '\t'.join(head)
        else:
            p = []
        return '\n'.join(p)

    def set_peek(self, dataset, **kwd):
        """
        expects a .pheno file in the extra_files_dir - ugh
        note that R is weird and does not include the row.name in
        the header. why?"""
        if not dataset.dataset.purged:
            pp = os.path.join(dataset.extra_files_path, '%s.pheno' % dataset.metadata.base_name)
            try:
                p = open(pp, 'r').readlines()
            except Exception:
                p = ['##failed to find %s' % pp, ]
            dataset.peek = ''.join(p[:5])
            dataset.blurb = 'Galaxy Rexpression composite file'
        else:
            dataset.peek = 'file does not exist\n'
            dataset.blurb = 'file purged from disk'

    def get_peek(self, dataset):
        """
        expects a .pheno file in the extra_files_dir - ugh
        """
        pp = os.path.join(dataset.extra_files_path, '%s.pheno' % dataset.metadata.base_name)
        try:
            p = open(pp, 'r').readlines()
        except Exception:
            p = ['##failed to find %s' % pp]
        return ''.join(p[:5])

    def get_file_peek(self, filename):
        """
        can't really peek at a filename - need the extra_files_path and such?
        """
        h = '## rexpression get_file_peek: no file found'
        try:
            h = open(filename, 'r').readlines()
        except Exception:
            pass
        return ''.join(h[:5])

    def regenerate_primary_file(self, dataset):
        """
        cannot do this until we are setting metadata
        """
        bn = dataset.metadata.base_name
        flist = os.listdir(dataset.extra_files_path)
        rval = ['<html><head><title>Files for Composite Dataset %s</title></head><p/>Comprises the following files:<p/><ul>' % (bn)]
        for i, fname in enumerate(flist):
            sfname = os.path.split(fname)[-1]
            rval.append('<li><a href="%s">%s</a>' % (sfname, sfname))
        rval.append('</ul></html>')
        with open(dataset.file_name, 'w') as f:
            f.write("\n".join(rval))
            f.write('\n')

    def init_meta(self, dataset, copy_from=None):
        if copy_from:
            dataset.metadata = copy_from.metadata

    def set_meta(self, dataset, **kwd):
        """
        NOTE we apply the tabular machinary to the phenodata extracted
        from a BioC eSet or affybatch.

        """
        Html.set_meta(self, dataset, **kwd)
        try:
            flist = os.listdir(dataset.extra_files_path)
        except Exception:
            if verbose:
                gal_Log.debug('@@@rexpression set_meta failed - no dataset?')
            return False
        bn = dataset.metadata.base_name
        if not bn:
            for f in flist:
                n = os.path.splitext(f)[0]
                bn = n
                dataset.metadata.base_name = bn
        if not bn:
            bn = '?'
            dataset.metadata.base_name = bn
        pn = '%s.pheno' % (bn)
        pp = os.path.join(dataset.extra_files_path, pn)
        dataset.metadata.pheno_path = pp
        try:
            pf = open(pp, 'r').readlines()  # read the basename.phenodata in the extra_files_path
        except Exception:
            pf = None
        if pf:
            h = pf[0].strip()
            h = h.split('\t')  # hope is header
            h = [escape(x) for x in h]
            dataset.metadata.column_names = h
            dataset.metadata.columns = len(h)
            dataset.peek = ''.join(pf[:5])
        else:
            dataset.metadata.column_names = []
            dataset.metadata.columns = 0
            dataset.peek = 'No pheno file found'
        if pf and len(pf) > 1:
            dataset.metadata.pheCols = self.get_phecols(phenolist=pf)
        else:
            dataset.metadata.pheCols = [('', 'No useable phenotypes found', False), ]
        if not dataset.info:
            dataset.info = 'Galaxy Expression datatype object'
        if not dataset.blurb:
            dataset.blurb = 'R loadable BioC expression object for the Rexpression Galaxy toolkit'
        return True

    def make_html_table(self, pp='nothing supplied from peek\n'):
        """
        Create HTML table, used for displaying peek
        """
        out = ['<table cellspacing="0" cellpadding="3">', ]
        try:
            # Generate column header
            p = pp.split('\n')
            for i, row in enumerate(p):
                lrow = row.strip().split('\t')
                if i == 0:
                    orow = ['<th>%s</th>' % escape(x) for x in lrow]
                    orow.insert(0, '<tr>')
                    orow.append('</tr>')
                else:
                    orow = ['<td>%s</td>' % escape(x) for x in lrow]
                    orow.insert(0, '<tr>')
                    orow.append('</tr>')
                out.append(''.join(orow))
            out.append('</table>')
            out = "\n".join(out)
        except Exception as exc:
            out = "Can't create html table %s" % str(exc)
        return out

    def display_peek(self, dataset):
        """
        Returns formatted html of peek
        """
        out = self.make_html_table(dataset.peek)
        return out


class Affybatch(RexpBase):
    """
    derived class for BioC data structures in Galaxy
    """

    file_ext = "affybatch"

    def __init__(self, **kwd):
        RexpBase.__init__(self, **kwd)
        self.add_composite_file('%s.affybatch',
                                description='AffyBatch R object saved to file',
                                substitute_name_with_metadata='base_name', is_binary=True)


class Eset(RexpBase):
    """
    derived class for BioC data structures in Galaxy
    """
    file_ext = "eset"

    def __init__(self, **kwd):
        RexpBase.__init__(self, **kwd)
        self.add_composite_file('%s.eset',
                                description='ESet R object saved to file',
                                substitute_name_with_metadata='base_name', is_binary=True)


class MAlist(RexpBase):
    """
    derived class for BioC data structures in Galaxy
    """
    file_ext = "malist"

    def __init__(self, **kwd):
        RexpBase.__init__(self, **kwd)
        self.add_composite_file('%s.malist',
                                description='MAlist R object saved to file',
                                substitute_name_with_metadata='base_name', is_binary=True)


class LinkageStudies(Text):
    """
    superclass for classical linkage analysis suites
    """
    test_files = [
        'allegro_descent', 'allegro_fparam', 'allegro_ihaplo',
        'alohomora_gts', 'alohomora_maf', 'alohomora_map', 'alohomora_ped',
        'ghm_haplo', 'ghm_lod',
        'linkage_datain', 'linkage_map', 'linkage_pedin',
        'merlin_chr', 'merlin_flow', 'merlin_lod',
        'simwalk_hef', 'swift_lod'
    ]

    def __init__(self, **kwd):
        Text.__init__(**kwd)
        self.lcount = 0
        self.max_lines = 2000
        # iterate whole file without errors
        self.eof_res = True

    @staticmethod
    def tokenizer(line, sep=None):
        """
        General purpose string tokenizer
        """
        if sep is None:
            return LinkageStudies.tokenizer(line)

        return line.splitlines()[0].split(sep)

    def eof_function(self):
        """
        Overridable end-of-file function
        """
        return self.eof_res

    @staticmethod
    def __is_binary_file(fstream):
        """
        Private binary file tester
        """
        result = False
        if '\x00' in fstream.readline(512):
            result = True

        fstream.seek(0)
        return result

    def __per_line_op(self, line):
        """
        Private per-line operation and line counter
        """
        self.lcount += 1
        if self.lcount > self.max_lines:
            return False

        return self.line_op(line)

    def line_op(self, line):
        """
        Overridable per line operation
        """
        return None

    def header_check(self, fio):
        """
        Overrideable post-binary file check function
        """
        return True

    def sniff(self, filename):
        """
        >>> from galaxy.datatypes.sniff import get_test_fname
        >>> cname = self.__class__.__name__
        >>> file_true = "linkstudies." + eval(cname)().file_ext
        >>> fname = get_test_fname(file_true)
        >>> eval(cname)().sniff(fname)
        True

        >>> result_array = []
        >>> for file_false in LinkageStudies.test_files:
        >>>     fname = get_test_fname("linkstudies." + file_false)
        >>>     if fname != file_true:
        >>>         res = eval(cname)().sniff(fname)
        >>>         result_array.push(res)
        >>>
        >>> set(result_array)
        {False}
        """
        with open(filename, "r") as fio:

            if LinkageStudies.__is_binary_file(fio):
                return False

            if not self.header_check(fio):
                return False

            for line in fio:
                line_res = self.__per_line_op(line)

                if line_res != None:
                    return line_res

            return self.eof_function()
        return False


class PreMakePed(LinkageStudies):
    """
    Common linkage pedin format
    Extended linkage pedigree file containing:
         pedigree_id, individ_id, fath_id, moth_id, gender, affectation, [genotypes]
    """
    file_ext = "linkage_pedin"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.num_colns = None
        self.max_lines = 1000

    def line_op(self, line):

        tokens = line.splitlines()[0].strip().split()

        if self.num_colns is None:
            self.num_colns = len(tokens)

        elif self.num_colns != len(tokens):
            return False

        try:
            if set([int(val) > 0 for val in tokens]) != {True}:
                return False

        except ValueError:
            return False

        return None


class Pedfile(PreMakePed):
    """
    Truncated PreMakePed file:
       - individuals specified on a single line
       - no genotype data
    """
    file_ext = "alohomora_ped"

    def __init__(self, **kwd):
        PreMakePed.__init__(**kwd)
        self.num_colns = 6


class GenotypeMatrix(LinkageStudies):
    """
    Sample matrix of genotypes
    - GTs as columns
    """
    file_ext = "alohomora_gts"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        # modern GT chipsets max at 5M
        self.max_lines = 5000000
        self.num_cols = -1

    def line_op(self, line):
        tokens = line.split('\t')

        if self.num_cols == -1:
            self.num_cols = len(tokens)

        elif self.num_cols != len(tokens):
            return False

        if not VALID_GENOTYPES_LINE.match(line):
            return False

        return None


    def header_check(self, fio):
        header_elems = fio.readline().splitlines()[0].strip().split('\t')

        try:
            [int(sid) > 0 for sid in header_elems[1:]]
        except ValueError:
            return False

        return True


class MarkerMap(LinkageStudies):
    """
    Map of genetic markers including physical and genetic distance
    Common input format for linkage programs

    chrom, genetic pos, markername, physical pos, Nr
    """
    file_ext = "linkage_map"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        # sensible linkage should not exceed more than 500,000 markers
        self.max_lines = 500000

    def header_check(self, fio):
        headers = fio.readline().split()

        if len(headers) == 5 and headers[0].lower()[:4] == "#Chr":
            return True

        return False

    def line_op(self, line):

        try:
            chrm, gpos, nam, bpos, row = LinkageStudies.tokenizer(line)

            float(gpos)
            int(bpos)

            try:
                int(chrm)
            except ValueError:
                if not chrm.lower()[0] in ('x', 'y', 'm'):
                    return False

        except ValueError:
            return False

        return None


class AlohomoraMarkerMap(LinkageStudies):
    """
    Map of genetic markers including physical and genetic distance

    chrom, markername, genetic pos, physical pos, markername, "x"
    """
    file_ext = "alohomora_map"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        # sensible linkage should not exceed more than 500,000 markers
        self.max_lines = 500000

    def header_check(self, fio):
        headers = fio.readline().split()

        if len(headers) == 6 and headers[0].lower()[:3] == "chr":
            return True

        return False

    def line_op(self, line):
        try:
            chrm, nam1, gpos, bpos, nam2, junk = LinkageStudies.tokenizer(line)

            float(gpos)
            int(bpos)

            try:
                int(chrm)
            except ValueError:
                if not chrm.lower()[0] in ('x', 'y', 'm'):
                    return False

            if nam1 != nam2:
                return False

        except ValueError:
            return False

        return None


class AlohomoraMAF(LinkageStudies):
    """
    Minor Allele Frequencies of marker lists
    """
    file_ext = "alohomora_maf"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.max_lines = 5000000

    def header_check(self, fio):
        fio.readline() # seek ahead
        return True

    def line_op(self, line):
        tokens = line.split()

        if len(tokens) != 2:
            return False

        try:
            int(tokens[0])
            return False
        except ValueError:
            pass

        try:
            float(tokens[1])
        except ValueError:
            return False

        return None


class DataIn(LinkageStudies):
    """
    Common linkage input file for intermarker distances
    and recombination rates
    """
    file_ext = "linkage_datain"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.num_markers = None
        self.intermarkers = 0

    def eof_function(self):
        return (self.num_markers - 1) == self.intermarkers

    def line_op(self, line):

        tokens = line.splitlines()[0].strip().split()

        try:

            if self.lcount == 1:

                self.num_markers = int(tokens[0])
                map(int, tokens[1:])

            elif self.lcount == 2:

                map(float, tokens)

                if len(tokens) != 4:
                    return False

            elif self.lcount == 3:

                map(int, tokens)
                last_token = int(tokens[-1])

                if self.num_markers is None:
                    return False

                if len(tokens) != last_token:
                    return False

                if self.num_markers == last_token:
                    return False

            elif int(tokens[0]) == 3 and int(tokens[1]) == 2:
                self.intermarkers += 1

            return None

        except ValueError:
            return False


class AllegroHaplo(PreMakePed):
    """
    Allegro output format for phased haplotypes
    """
    file_ext = "allegro_ihaplo"

    def header_check(self, fio):
        header = []
        while True:
            line = fio.readline()
            if line.startswith("          "):
                header.append(line.splitlines()[0])
            else:
                break

        if header == []:
            return False

        # transpose headers
        markers = ["".join(x[::-1]).strip() for x in zip(*header)]
        markers = [mark for mark in markers if mark != ""]

        for mark in markers:
            if len(mark.split(" ")) > 0:
                return False


class AllegroDescent(AllegroHaplo):
    """
    Allegro output format for founder allele groups
    """
    file_ext = "allegro_descent"


class AllegroLOD(LinkageStudies):
    """
    Allegro output format for LOD scores
    """
    file_ext = "allegro_fparam"

    def header_check(self, fio):
        header = fio.readline().splitlines()[0].split()
        if not(header[0] == "family"
               and header[1] == "location"
               and header[2] == "LOD"
               and header[3] == "marker"):
            return False
        return True

    def line_op(self, line):
        tokens = LinkageStudies.tokenizer(line)

        try:
            int(tokens[0])
            float(tokens[1])

            if tokens[2] != "-inf":
                float(tokens[2])

        except ValueError:
            return False

        return None


class AllegroSetup(LinkageStudies):
    """
    Allegro input setup file
    """
    file_ext = "allegro_in"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.max_lines = 100
        self.find_line = {
            'PREFILE' : False,
            'DATFILE' : False,
            'MODEL' : False
        }
        self.eof_res = False

    def line_op(self, line):
        for f_line in self.find_line:
            if line.startswith(f_line):
                self.find_line[f_line] = True

        if set(self.find_line.values()) == {True}:
            return True

        return None


class GHMHaplo(LinkageStudies):
    """
    Genehunter output phased haplotypes
    """
    file_ext = "ghm_haplo"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.num_markers = -1

    def header_check(self, fio):
        return fio.readline().startswith("*****")


    def line_op(self, line):
        tokens = LinkageStudies.tokenizer(line)

        try:
            num_tokens = len(map(int, tokens))

            if self.num_markers == -1:
                self.num_markers = num_tokens - 4

            elif self.num_markers != num_tokens:
                self.num_markers = -1
                return False

        except ValueError:
            return False

        return None


class GHMLOD(LinkageStudies):
    """
    Genehunter output file for parametric LOD
    """
    file_ext = "ghm_lod"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.lod_header = -1
        self.eof_res = False
        self.lin_fin = "position  LOD score    NPL score  p-value    information"

    def line_op(self, line):

        if line.startswith(self.lin_fin):
            self.lod_header = self.lcount
            return None

        if self.lod_header != -1:
            if self.lod_header < self.lcount <= self.lod_header + 20:
                tokens = LinkageStudies.tokenizer(line)

                try:
                    map(float, tokens)
                    if self.lod_header + 20:
                        return True

                except ValueError:
                    return False

        return None


class MerlinLOD(GHMLOD):
    """
    Merlin output file for parametric LOD
    """
    file_ext = "merlin_lod"

    def __init__(self, **kwd):
        GHMLOD.__init__(**kwd)
        self.lin_fin = "       POSITION        LOD      ALPHA       HLOD"

    def header_check(self, fio):
        return fio.readline().split()[0] == "MERLIN"


class MerlinChr(LinkageStudies):
    """
    Merlin output file for phased haplotypes
    """
    file_ext = "merlin_chr"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.indiv_line = -1

    def header_check(self, fio):
        return fio.readline().split()[0] == "FAMILY"

    def __internal_line_op(self, alleles):
        """
        Merlin-specific data handler for both alleles and
        later founder groups
        """
        try:
            map(int, alleles)

            if self.indiv_line + 20:
                return True

        except ValueError:
            return False

        return None


    def line_op(self, line):

        if line.contains("(F)") or line.contains("(M)"):
            self.indiv_line = self.lcount
            return None

        if self.indiv_line != -1:
            if self.indiv_line < self.lcount <= self.indiv_line + 20:
                tokens = LinkageStudies.tokenizer(line)
                alleles = filter(
                    [alle for alle in tokens if alle not in ('?', ':', ',', '|')]
                )

                return self.__internal_line_op(alleles)

        return None


class MerlinHaplo(MerlinChr):
    """
    Merlin output file for descent data (founder allele groups)
    """
    file_ext = "merlin_flow"

    def __internal_line_op(self, alleles):
        try:
            if set(map([alle for alle in alleles if ord(alle) in range(65, 90)])) != {True}:
                return False

        except ValueError:
            return False

        return None


class SwiftLOD(LinkageStudies):
    """
    Swiftlink output LOD file
    """
    file_ext = "swift_lod"

    def header_check(self, fio):
        return fio.readline().startswith("marker  position        lod")

    def line_op(self, line):

        tokens = LinkageStudies.tokenizer(line)

        if tokens[0] == "-":
            if len(tokens) != 3:
                return False

            try:
                float(tokens[1])
                float(tokens[2])
            except ValueError:
                return False

            return None

        # marker name
        if len(tokens) != 2:
            return False

        try:
            float(tokens[1])
        except ValueError:
            return False

        return None


class SimwalkHEF(LinkageStudies):
    """
    Simwalk output HEF file. Contains LOD, phased haplotypes,
    and founder groups
    """
    file_ext = "simwalk_hef"

    def __init__(self, **kwd):
        LinkageStudies.__init__(**kwd)
        self.lod_header = -1
        self.hap_header = -1

    def header_check(self, fio):
        return fio.readline().split()[0] == "HEF"

    def eof_function(self):
        return self.lod_header != -1 or self.hap_header != -1

    def line_op(self, line):

        if line.startswith("Marker  Map: Haldane-cM  Number of  Allele Names and"):
            self.hap_header = -1
            self.lod_header = self.lcount

        elif line.startswith("Num. of Individuals                             Pheno"):
            self.hap_header = self.lcount
            self.lod_header = -1

        try:
            if self.lod_header != -1:
                if self.lod_header + 2 < self.lcount < self.lod_header + 12:

                    tokens = LinkageStudies.tokenizer(line)

                    if (self.lcount - self.lod_header) % 2 == 0:
                        # even lines, markers
                        if line.startswith("          "):
                            return False

                        float(tokens[1])
                        float(tokens[2])
                        int(tokens[3])

                    else:
                        # odd lines, intermarker distances
                        if not line.startswith("          "):
                            return False

                        int(tokens[0])
                        float(tokens[1])
                        int(tokens[2])
                        float(tokens[3])


            elif self.hap_header != -1:
                tokens = LinkageStudies.tokenizer(line)
                map(int, tokens)

                # linter prefers less nesting...
                if self.lcount == self.hap_header + 6 and len(tokens) != 5:
                    return False

                elif ((self.hap_header + 6 < self.lcount < self.hap_header + 16)
                      and len(tokens) != 6):
                    return False

        except ValueError:
            return False

        return None


if __name__ == '__main__':
    import doctest
    doctest.testmod(sys.modules[__name__])
