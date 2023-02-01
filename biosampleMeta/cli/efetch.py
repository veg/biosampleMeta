import argparse

from biosampleMeta.scrape import efetch


def efetch_cli():
    parser = argparse.ArgumentParser(description='Perform an NCBI efetch.')
    parser.add_argument('--db', type=str, help='database to access')
    parser.add_argument('--accession', type=str, help='accessions to download')
    parser.add_argument('--email', type=str, help='email address')
    args = parser.parse_args()
    Entrez.email = args.email.replace('__at__', '@')
    with open(args.accession) as accession_file:
        accessions = [a.strip() for a in accession_file.readlines()]
    for accession in accessions:
        sleep(1)
        output_json = 'output/%s.json' % accession
        efetch_serialize(args.db, accession, 'output/temp.xml', output_json)



