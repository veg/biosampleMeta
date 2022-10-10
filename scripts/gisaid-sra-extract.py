import csv
import os
import re

from bs4 import BeautifulSoup


base_dir = 'output/sra-parsing'
xmls = [
    f
    for f in os.listdir(base_dir)
    if f.split('.')[-1] == 'xml'
]

gisaid_accession_regex = re.compile('EPI_ISL_\d+')

csv_file = open('output/sra-gisaid.csv', 'w')
writer = csv.writer(csv_file)
writer.writerow(['sra_accession', 'gisaid_accession'])
for xml_filename in xmls:
    with open('%s/%s' % (base_dir, xml_filename)) as xml_file:
        xml = xml_file.read()
        soup = BeautifulSoup(xml, 'xml')
        try:
            result = soup.find('RUN')
            sra_accession = result.get('accession')
            gisaid_accession = gisaid_accession_regex.search(xml).group()
            writer.writerow([sra_accession, gisaid_accession])
            print('wrote ', xml_filename)
        except:
            pass
csv_file.close()
