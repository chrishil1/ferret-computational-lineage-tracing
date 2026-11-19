# Sample sheet
The sample sheet should include information about all of your samples in a tsv format.

Column 1 - Sample name
Column 2 - comma separated list of R1 fastq files, the R2 must have the same name execpt R2 instead of R1
Column 3 - Sample wells.
	Wells are specified in blocks, ranges, or individually like this:
        'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.
        'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.
        'C4' specifies a single well.
        Multiple selections are joined by commas (NO SPACES), e.g. 'A1-A6,B1:D3,C4'
Column 4 - Path to the raw data. If multiple fastq files are supplied and are in different directories, you may pass a comma separated list of directories the same length and order as the fastq files.
Column 5 - The names associated with the seprate fastq files, example run1, run2. All unique fastq files must have a unique name and all of the fastq files that are the same must have the same name.