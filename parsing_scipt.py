#Instrantiate dictionary to hold reads
mate_pairs = {}


# The list of query names we wish to select for.
QNAME_list = ["SRR1972739.1", "SRR1972739.2", "SRR1972739.3", "SRR1972739.4", "SRR1972739.5", "SRR1972739.6", "SRR1972739.7", "SRR1972739.8", "SRR1972739.9", "SRR1972739.10"]


# make the checking a bit faster by making the list a set
QNAME_list = set(QNAME_list) 


# Open and read the SAM file created previously
# Read each line and store the queryname of each line in QNAME
# If QNAME is in QNAME_list, store SAM file alignment in mate_pairs dictionary
with open("bwa.sam", 'r') as samFile:
    for line in samFile:
        QNAME = line.split('\t')[0]
        if QNAME in QNAME_list:
            try:
                mate_pairs[QNAME].append(line)
            except KeyError:
                mate_pairs[QNAME] = [line]


# Open and write to new parsed SAM file.
# Iterate through the items of the mate_pairs dictionary.
# Write the values of the dictionary (i.e. the alignments) to the new parsed SAM file.
with open("parsed_SAM_file.sam", "w") as outfile:
    for read,mates in mate_pairs.items():
        for mate in mates:
            outfile.write('%s\n' % (mate))
