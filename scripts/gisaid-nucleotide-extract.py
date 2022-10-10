import csv
import os
import re

from bs4 import BeautifulSoup


base_dir = 'output/nucleotide-parsing'
xmls = [
    f
    for f in os.listdir(base_dir)
    if f.split('.')[-1] == 'xml'
]

gisaid_accession_regex = re.compile('EPI_ISL_\d+')

csv_file = open('output/nucleotide-gisaid.csv', 'w')
writer = csv.writer(csv_file)
writer.writerow(['nucleotide_accession', 'gisaid_accession'])
for xml_filename in xmls:
    with open('%s/%s' % (base_dir, xml_filename)) as xml_file:
        xml = xml_file.read()
        soup = BeautifulSoup(xml, 'xml')
        result = soup.find('Textseq-id_accession')
        try:
            nucleotide_accession = result.text
            gisaid_accession = gisaid_accession_regex.search(xml).group()
            writer.writerow([nucleotide_accession, gisaid_accession])
        except:
            print('issue with', xml_filename)
csv_file.close()
