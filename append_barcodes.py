from getopt import getopt
import sys
import os
from glob import glob
import gzip
from types import StringTypes

HELP_STRING = """
Given a read file (or paired read files) and a barcode file, appends the barcode
to the read name in a manner compatible with nudup processing.
    -1   Forward read file
    -2   Reverse read file
    -b   Barcode (UMI) file
If a -2 parameter is not provided then it will be assumed that the data is
single ended. The new files will be named for the original files + "_bc.fq".
VERSION: $Revision: 0.01 $
"""

def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        if (fh.endswith(".gz")):
            fh = gzip.open(fh,'r')
        else:
            fh = file(fh)

    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise Exception("Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry.")
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )

        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))

def appendBarcode(fwdTitle, fwdSeq, fwdQual,
                  revTitle, revSeq, revQual,
                  barcode, bcLength, outFwd, outRev):
    bc = barcode[-bcLength:]
    sys.stderr.write(bc)
    outFwd.write("@%s:%s\n%s\n+\n%s\n" % (fwdTitle, bc, fwdSeq, fwdQual))
    if revSeq != None: # if this is a paired end read
        outRev.write("@%s:%s\n%s\n+\n%s\n" % (revTitle, bc, revSeq, revQual))

# retrieve the user parameters
fwdFilename = None
revFilename = None
barcodeFilename = None
bcLength = 6

try:
    optlist, args = getopt(sys.argv[1:], "h1:2:b:l:")
except:
    print "Error retrieving options"
    print ""
    print HELP_STRING
    sys.exit(1)

for (opt, opt_arg) in optlist:
    if opt == "-h":
        print ""
        print HELP_STRING
        sys.exit(1)
    elif opt == "-1":
        fwdFilename = opt_arg
    elif opt == "-2":
        revFilename = opt_arg
    elif opt == "-b":
        barcodeFilename = opt_arg
    elif opt == "-l`":
        bcLength = opt_arg

if fwdFilename == None:
    print "\nYou must provide a fwd FASTQ filename or both fwd & rev filenames."
    print
    print HELP_STRING
    sys.exit(1)

if barcodeFilename == None:
    print "\nYou must provide a barcode FASTQ filename."
    print
    print HELP_STRING
    sys.exit(1)

recCount = 1

fwdIt = FastqIterator(fwdFilename)
(fwdTitle, fwdSeq, fwdQual) = fwdIt.next()
fwdReadName = fwdTitle.split()[0]

bcIt = FastqIterator(barcodeFilename)
(bcTitle, bcSeq, bcQual) = bcIt.next()
bcReadName = bcTitle.split()[0]

if revFilename:
    revIt = FastqIterator(revFilename)
    (revTitle, revSeq, revQual) = revIt.next()
    revRoot = os.path.splitext(revFilename)[0]
    revReadName = revTitle.split()[0]
    if (revFilename.endswith(".gz")):
        outRev = gzip.open(revRoot + "_trimmed.fq.gz", "wb")
    else:
        outRev = open(revRoot + "_trimmed.fq", "w")
else:
    revIt = None
    revReadName = revTitle = revSeq = revQual = None
    outRev = None

if fwdReadName != bcReadName:
    print "The fwd read name and barcode read name don't match"
    print "fwd name: '%s'" % fwdReadName
    print "bc name: '%s'" % bcReadName
    raise Exception("fwd and barcode names don't match.")

if revReadName and fwdReadName != revReadName:
    print "The fwd read name and rev read name don't match"
    print "fwd name: '%s'" % fwdReadName
    print "rev name: '%s'" % revReadName
    raise Exception("fwd and rev names don't match.  Not paired end sequence.")

if len(bcSeq) < bcLength:
    print "The barcode length doesn't match the actual barcode sequence"
    print "Barcode sequence: '%s'" % (bcSeq.length)
    print "Barcode length  : '%d'" % bcLength
    raise Exception("The barcode length doesn't match the actual barcode sequence")

fwdRoot = os.path.splitext(fwdFilename)[0]
if (fwdFilename.endswith(".gz")):
    outFwd = gzip.open(fwdRoot + "_bc.fq.gz", "wb")
else:
    outFwd = open(fwdRoot + "_bc.fq", "w")

while True:
    # trim one record and add to output
    appendBarcode(fwdTitle, fwdSeq, fwdQual, revTitle,
        revSeq, revQual, bcSeq, bcLength, outFwd, outRev)

    try:
        (fwdTitle, fwdSeq, fwdQual) = fwdIt.next()
        (bcTitle, bcSeq, bcQual) = bcIt.next()
        if revFilename:
            (revTitle, revSeq, revQual) = revIt.next()
    except StopIteration:
        break

    recCount += 1

print "\tDone with this set. there were %s records" % (recCount)
