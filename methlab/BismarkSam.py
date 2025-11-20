class BismarkSam(object):
    """
    Class to import and manipulate the XM tags of SAM files exported from Bismark.
    """
    def __init__(self, read) -> None:
        """
        Initialise class function.
        """
        self.read = read
        if read[13][:5] != "XM:Z:":
            raise ValueError("Element 13 of this read is not an XM tag. This read does not appear to be a Bismark-sorted SAM file.")
        self.id = read[0]
        self.chr = read[2]
        self.length = len(read[9])
        self.seq = read[9]
        self.flag = read[1]
        self.xm_tag = self.get_tag("XM:Z:")
        self.strand = self.get_tag("YS:Z:")
    
    def get_tag(self, tag_name):
        """
        Retrieve a tag from the SAM entry with a specific tag.

        This checks each tag in the SAM entry for whether the input tag name
        matches the *start* of each tag, and pulls out the string after the end
        of the tag name.

        Parameters
        ----------
        tag_name : str
            Tag prefix to search for, for example ``"YS:Z:"`` or ``"XM:Z:"``.

        Returns
        -------
        str
            The value of the matching tag, or ``"NA"`` if no such tag is present.

        Raises
        ------
        ValueError
            If more than one tag is found with the same tag ID.

        Examples
        --------
        Strand can be coded as ``"YS:Z:CTOT"``, and you want to retrieve
        ``"CTOT"``:

        >>> read.get_tag("YS:Z:")
        'CTOT'
        """
        tag_length = len(tag_name)
        tag_value = [tag[tag_length:] for tag in self.read if tag[:tag_length] == tag_name]
        number_of_values = len(tag_value)
        if number_of_values == 0:
            return 'NA'
        elif number_of_values == 1:
            return tag_value[0]
        elif number_of_values > 1:
            raise ValueError("Multiple tags were found with the same tag ID")

    def count_mC(self):
        """
        Count the number of methylated and unmethylated cytosines on this read.

        Returns
        -------
        list of int
            Two-element list ``[methylated, unmethylated]``.
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
        Check whether all methylated cytosines appear together in this read.

        A cluster is defined as one or more consecutive methylated cytosines
        with no unmethylated cytosines in between. It is important that
        non-cytosines have been removed from the read before calling this.

        Returns
        -------
        bool
            ``True`` if there is at most one methylated cluster, ``False``
            if multiple clusters are found.
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
        Summarise methylation status for this read.

        Returns
        -------
        list
            A 5-element list ``[chromosome, methylated, unmethylated, length, cluster]``,
            where ``cluster`` is ``True``/``False`` or ``"NA"`` if not applicable.
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
    Import a SAM file and parse each read into a :class:`BismarkSam` object.

    Parameters
    ----------
    path : str
        Path to a SAM file exported from Bismark.

    Returns
    -------
    list of BismarkSam
        Parsed SAM reads wrapped in :class:`BismarkSam` objects.
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