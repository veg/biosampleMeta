import csv


csv_file = open('output/sra-gisaid.csv')
reader = csv.DictReader(csv_file)
gisaid2sra = {
    row['gisaid_accession']: row['sra_accession']
    for row in reader
}
with open('input/gisaid-intersection.txt') as f:
    sra_accessions = [gisaid2sra[l.strip()] for l in f.readlines()]

with open('output/sra-accessions-intersecting-with-gisaid.txt', 'w') as f:
    f.write('\n'.join(sra_accessions))
