class BismarkSam(object):
    """
    Class to import and manipulate the XM tags of SAM files exported from Bismark.
    """
    def __init__(self, read) -> None:
        """
        Initialise class function
        """
        if read[13][:5] != "XM:Z:":
            raise ValueError("Element 13 of this read is not an XM tag. This read does not appear to be a Bismark-sorted SAM file.")
        self.id = read[0]
        self.chr = read[2]
        self.length = len(read[9])
        self.xm_tag = read[13][5:]
        self.seq = read[9]
    
    def count_mC(self):
        """
        Count how many cytosines on a read are methylated, unmethylated, and the total.
        """
        upper = 0
        lower = 0
        up="HXZ"
        lo="hxz"
        for i in self.xm_tag:
            if i in up:
                upper+=1
            elif i in lo:
                lower+=1
        return [upper, lower]

    def mC_cluster(self) :
        """
        Check whether all the methylated cytosines appear together in a read
        (there  are no unmethylated cytosines between methylated cytosines).
        
        It is important that non-cytosines have been removed from the read.
        """
        flag = False
        index = 0
        trimmed_tag = [c for c in self.xm_tag if (c.islower() or c.isupper()) ]
        n = len(trimmed_tag)
        
        # Check for clusters in reads with two or more cytosines.
        while index < n:
            if trimmed_tag[index].isupper():
                if (flag == True) :
                    return False
                while index < n and trimmed_tag[index].isupper():
                    index += 1
                flag = True
            else :
                index += 1
        return True
    
    def mC_per_read(self):
        """
        Summarise the number of methylated reads, unmethylated reads, total read
        length, and whether or not cytosines are occur next to one another
        """
        total = self.length
                        
        mC, uC = self.count_mC()
        if mC > 1:
            cluster = self.mC_cluster()
        else:
            cluster = 'NA'
        
        return [self.chr, mC, uC, total, cluster]

def read_SAM(input:str):
    """
    Import a SAM file.

    Parameters
    ==========
    input: str
        Path to a SAM file.
    
    Returns
    =======
    A list of BismarkSam objects.
    """
    # Import samfile line by line
    reads = [line.split("\t") for line in open(input, "r").read().split("\n")]
    print(input)
    print("File has {} lines".format( len(reads)) )
    # Pull out header information
    header = [read for read in reads if read[0].startswith('@')]
    header = ['\t'.join(read) + '\n' for read in header]
    # Keep reads if they are not empty and do not start with '@'
    reads = [read for read in reads if (read != ['']) & ( read[0].startswith("@") is False ) ]

    reads = [BismarkSam(read) for read in reads]

    return reads