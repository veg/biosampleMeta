import os
import re
from time import sleep

from Bio import Entrez
from bs4 import BeautifulSoup


Entrez.email = os.environ.get('ENTREZ_EMAIL') or None
Entrez.api_key = os.environ.get('ENTREZ_API_KEY') or None


wildcard_constraints:
  bp_accession="[^/]+"

with open('input/tax_ids.txt') as f:
  tax_ids = f.read().splitlines()

def efetch(db, accession):
  handle = Entrez.efetch(db=db, id=accession)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def elink(db, dbfrom, id_, retmax=1000):
  handle = Entrez.elink(db=db, dbfrom=dbfrom, id=id_, retmax=retmax)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def esearch(db, term, retmax=1000):
  handle = Entrez.esearch(db=db, term=term, retmax=retmax)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def read_soup(filename):
    with open(filename) as f:
      soup = BeautifulSoup(f, 'xml')
    return soup


def write_soup(soup, filename):
  pretty = soup.prettify()
  with open(filename, 'w') as f:
      f.write(pretty)


def write_ids(ids, filename):
  with open(filename, 'w') as f:
    f.write('\n'.join(ids))


rule bioproject_ids_from_tax_id:
  output:
    xml="output/{tax_id}/search.xml",
    accessions="output/{tax_id}/bioproject_accessions.txt"
  run:
    tax_id = 'txid' + wildcards.tax_id
    soup = esearch('bioproject', tax_id)
    write_soup(soup, output.xml)

    ids = [
      id_.text.strip()
      for id_ in soup.find('IdList').findAll('Id')
    ]
    write_ids(ids, output.accessions)

rule all_bioproject_ids:
  input:
    expand(
      "output/{tax_id}/bioproject_accessions.txt",
      tax_id=tax_ids
    )

rule bioproject_xml_fetch:
  output:
    "output/{tax_id}/{bp_accession}/bioproject.xml"
  run:
    soup = efetch('bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_query:
  input:
    rules.bioproject_xml_fetch.output[0]
  output:
    "output/{tax_id}/{bp_accession}/sra_accessions.xml"
  run:
    soup = elink('sra', 'bioproject', id_=wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_ids:
  input:
    rules.bioproject_sra_query.output[0]
  output:
    "output/{tax_id}/{bp_accession}/sra_accessions.txt"
  run:
    soup = read_soup(input[0])
    ids = [
      id_.text
      for id_ in soup.find('LinkSetDb').findAll('Id')
    ]
    write_ids(ids, output.accessions)

rule bioproject_biosample_query:
  input:
    rules.bioproject_xml_fetch.output[0]
  output:
    "output/{tax_id}/{bp_accession}/biosample_links.xml"
  run:
    soup = elink('biosample', 'bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_ids:
  input:
    rules.bioproject_biosample_query.output[0]
  output:
    "output/{tax_id}/{bp_accession}/biosample_accessions.txt"
  run:
    soup = read_soup(input[0])
    samn_ids = [
      id_.text.strip()
      for id_ in soup.find('LinkSetDb').findAll('Id')
    ]
    write_ids(samn_ids, output[0])
 
rule bioproject_biosample_entry:
  output:
    "output/{tax_id}/{bp_accession}/biosamples/{bs_accession}/entry.xml"
  run:
    soup = efetch('biosample', wildcards.bs_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_geo_accession:
  input:
    rules.bioproject_biosample_entry.output[0]
  output:
    "output/{tax_id}/{bp_accession}/biosamples/{bs_accession}/geo_accession.txt"
  run:
    soup = read_soup(input[0])
    geo_accession = soup.find('Id', {'db': "GEO"}).text.strip()
    write_ids([geo_accession], output[0])

rule bioproject_biosample_geo_xml:
  input:
    rules.bioproject_biosample_geo_accession.output[0]
  output:
    "output/{tax_id}/{bp_accession}/biosamples/{bs_accession}/geo.xml"
  shell:
    """
      sleep 3;
      GEO_ACCESSION=$(cat {input})
      wget "https://www.ncbi.nlm.nih.gov/geo/tools/geometa.cgi?acc=$GEO_ACCESSION&scope=full&mode=miniml" -O {output}
    """

#def get_bs_ids(wildcards):
#  filename = "bioprojects/%s/biosample_links.txt" % wildcards.bp_accession
#  with open(filename) as f:
#    bs_ids = [line.strip() for line in f.readlines()]
#  bs_files = [
#    "bioprojects/%s/biosamples/%s/entry.xml" % (wildcards.bp_accession, bs_id)
#    for bs_id in bs_ids
#  ]
#  return bs_files 
#
#
#rule all_bioproject_biosample_entry:
#  input:
#    "bioprojects/{bp_accession}/biosample_links.txt",
#    get_bs_ids
#  output:
#    'bioprojects/{bp_accession}/biosamples/all.txt'
#  run:
#    with open(output[0], 'w') as f:
#      f.write('success')
#
#rule all_biosample_links:
#  input:
#    expand(
#      'bioprojects/{bp_accession}/biosample_links.txt',
#      bp_accession=bp_accessions
#    )
#
#rule all_biosamples:
#  input:
#    expand(
#      'bioprojects/{bp_accession}/biosamples/all.txt',
#      bp_accession=ran_bp_accessions
#    )
