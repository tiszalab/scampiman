import subprocess
from subprocess import Popen, PIPE, STDOUT

subprocess.Popen(['which', 'samtools'])

subprocess.Popen([
    'samtools', 
    'ampliconstats', 
    '-o', 
    'scampi_test1/testtest.ampliconstats.tsv', 
    'SARS-CoV-2.primer.bed', 
    '/var/tmp/62818/scampi2/ts1_bar1_1a.ampclip.bam'
    ]
    )