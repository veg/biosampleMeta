import os
import re
from time import sleep

from Bio import Entrez
from bs4 import BeautifulSoup


Entrez.email = os.environ.get('ENTREZ_EMAIL') or None
Entrez.api_key = os.environ.get('ENTREZ_API_KEY') or None


wildcard_constraints:
  bp_accession="[^/]+"


def efetch(db, accession):
  handle = Entrez.efetch(db=db, id=accession)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def elink(db, dbfrom, accession):
  accession_id = re.sub('\D', '', accession)
  handle = Entrez.elink(db=db, dbfrom=dbfrom, id=accession_id, retmax=100)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def esearch(db, term, retmax=10):
  handle = Entrez.esearch(db=db, term=term, retmax=retmax)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def write_soup(soup, filename):
  pretty = soup.prettify()
  with open(filename, 'w') as f:
      f.write(pretty)


rule bioproject_xml_fetch:
  output:
    "bioprojects/{bp_accession}/entry.xml"
  run:
    soup = efetch('bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

#rule bioproject_sra_accessions:
#  input:
#    rules.bioproject_xml_fetch.output[0]
#  output:
#    "bioprojects/{bp_accession}/sra_accessions.txt"
#  run:

rule bioproject_biosample_link:
  output:
    "bioprojects/{bp_accession}/biosample_links.xml"
  run:
    soup = elink('biosample', 'bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_ids:
  input:
    rules.bioproject_biosample_link.output[0]
  output:
    "bioprojects/{bp_accession}/biosample_links.txt"
  run:
    with open(input[0]) as f:
      soup = BeautifulSoup(f, 'xml')
    samn_ids = [id_.text.strip() for id_ in soup.find('LinkSetDb').findAll('Id')]
    with open(output[0], 'w') as f:
      f.write('\n'.join(samn_ids))
 
rule bioproject_biosample_entry:
  output:
    "bioprojects/{bp_accession}/biosamples/{bs_accession}/entry.xml"
  run:
    soup = efetch('biosample', wildcards.bs_accession)
    write_soup(soup, output[0])

def get_bs_ids(wildcards):
  filename = "bioprojects/%s/biosample_links.txt" % wildcards.bp_accession
  with open(filename) as f:
    bs_ids = [line.strip() for line in f.readlines()]
  bs_files = [
    "bioprojects/%s/biosamples/%s/entry.xml" % (wildcards.bp_accession, bs_id)
    for bs_id in bs_ids
  ]
  return bs_files 


rule all_bioproject_biosample_entry:
  input:
    "bioprojects/{bp_accession}/biosample_links.txt",
    get_bs_ids
  output:
    'bioprojects/{bp_accession}/biosamples/all.txt'
  run:
    with open(output[0], 'w') as f:
      f.write('success')

rule bioproject_biosample_geo_accession:
  input:
    rules.bioproject_biosample_entry.output[0]
  output:
    "bioprojects/{bp_accession}/biosamples/{bs_accession}/geo_accession.txt"
  run:
    with open(input[0]) as f:
      soup = BeautifulSoup(f, 'xml')
    geo_accession = soup.find('Id', {'db': "GEO"}).text.strip()
    with open(output[0], 'w') as f:
      f.write(geo_accession)

rule bioproject_biosample_geo_xml:
  input:
    rules.bioproject_biosample_geo_accession.output[0]
  output:
    "bioprojects/{bp_accession}/biosamples/{bs_accession}/geo.xml"
  shell:
    """
      sleep 3;
      GEO_ACCESSION=$(cat {input})
      wget "https://www.ncbi.nlm.nih.gov/geo/tools/geometa.cgi?acc=$GEO_ACCESSION&scope=full&mode=miniml" -O {output}
    """
