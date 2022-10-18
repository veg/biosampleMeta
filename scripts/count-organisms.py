import os
from collections import Counter

from bs4 import BeautifulSoup


with open('output/sra-accessions-intersecting-with-gisaid.txt') as f:
    sra_accessions = set([l.strip() for l in f.readlines()])
tax_counter = Counter()
base_dir = 'output/sra-parsing'
xml_filenames = [fn for fn in os.listdir(base_dir) if fn.split('.')[-1] == 'xml']
for xml_filename in xml_filenames:
    xml_filepath = '%s/%s' % (base_dir, xml_filename)
    with open(xml_filepath) as xml_file:
        soup = BeautifulSoup(xml_file, 'xml')
    sra_accession = soup.find('RUN').get('accession')
    if sra_accession in sra_accessions:
        tax_counter[soup.find('TAXON_ID').text.strip()] += 1
print(tax_counter)
