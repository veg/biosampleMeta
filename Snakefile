import os
import re
import csv
import json
from time import sleep

from Bio import Entrez
from bs4 import BeautifulSoup
import xmltodict


Entrez.email = os.environ.get('ENTREZ_EMAIL') or None
Entrez.api_key = os.environ.get('ENTREZ_API_KEY') or None


wildcard_constraints:
  bp_accession="[^/]+"


lineage = {
  11269: 'Marburg marburgvirus|Marburgvirus|Filoviridae|Mononegavirales|Monjiviricetes|Haploviricotina|Negarnaviricota|Orthornavirae|Riboviria|Viruses',
  11320: 'Influenza A virus|Alphainfluenzavirus|Orthomyxoviridae|Articulavirales|Insthoviricetes|Polyploviricotina|Negarnaviricota|Orthornavirae|Riboviria|Viruses'
}


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


def read_ids(filename):
  with open(filename) as f:
    ids = f.read().splitlines()
  return ids


def write_ids(ids, filename):
  with open(filename, 'w') as f:
    f.write('\n'.join(ids))


def read_json(filename):
  with open(filename) as json_file:
    result = json.load(json_file)
  return result


def write_json(result, filename):
  with open(filename, 'w') as json_file:
    json.dump(result, json_file, indent=2)


def xml2json(xml_filename, json_filename):
  with open(xml_filename) as xml_file:
    xml = xml_file.read()
  converted = xmltodict.parse(xml)
  with open(json_filename, 'w') as json_file:
    json.dump(converted, json_file, indent=2)


rule bioproject_query:
  output:
    "output/{tax_id}/bioprojects/search.xml",
  run:
    tax_id = 'txid' + wildcards.tax_id
    soup = esearch('bioproject', tax_id)
    write_soup(soup, output[0])

rule bioproject_ids_from_query:
  input:
    rules.bioproject_query.output[0]
  output:
    "output/{tax_id}/bioprojects/ids.txt"
  run:
    soup = read_soup(input[0])
    ids = [
      id_.text.strip()
      for id_ in soup.find('IdList').findAll('Id')
    ]
    write_ids(ids, output[0])

rule all_bioproject_ids:
  input:
    expand(
      "output/{tax_id}/bioprojects/ids.txt",
      tax_id=read_ids('input/tax_ids.txt')
    )

rule bioproject_xml_fetch:
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/bioproject.xml"
  run:
    soup = efetch('bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_query:
  input:
    rules.bioproject_xml_fetch.output[0]
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/sra_accessions.xml"
  run:
    soup = elink('sra', 'bioproject', id_=wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_ids:
  input:
    rules.bioproject_sra_query.output[0]
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/sra_accessions.txt"
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
    "output/{tax_id}/bioprojects/{bp_accession}/biosample_links.xml"
  run:
    soup = elink('biosample', 'bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_ids:
  input:
    rules.bioproject_biosample_query.output[0]
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/biosample_accessions.txt"
  run:
    soup = read_soup(input[0])
    samn_ids = [
      id_.text.strip()
      for id_ in soup.find('LinkSetDb').findAll('Id')
    ]
    write_ids(samn_ids, output[0])
 
rule bioproject_biosample_entry:
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/entry.xml"
  run:
    soup = efetch('biosample', wildcards.bs_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_geo_accession:
  input:
    rules.bioproject_biosample_entry.output[0]
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/geo_accession.txt"
  run:
    soup = read_soup(input[0])
    geo_accession = soup.find('Id', {'db': "GEO"}).text.strip()
    write_ids([geo_accession], output[0])

rule bioproject_biosample_geo_xml:
  input:
    rules.bioproject_biosample_geo_accession.output[0]
  output:
    "output/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/geo.xml"
  shell:
    """
      sleep 3;
      GEO_ACCESSION=$(cat {input})
      wget "https://www.ncbi.nlm.nih.gov/geo/tools/geometa.cgi?acc=$GEO_ACCESSION&scope=full&mode=miniml" -O {output}
    """

rule assembly_query:
  output:
    "output/{tax_id}/assembly/search.xml"
  run:
    tax_id = 'txid' + wildcards.tax_id
    soup = esearch('assembly', tax_id)
    write_soup(soup, output[0])

rule assembly_ids_from_query:
  input:
    rules.assembly_query.output[0]
  output:
    "output/{tax_id}/assembly/ids.txt"
  run:
    soup = read_soup(input[0])
    ids = [
      id_.text.strip()
      for id_ in soup.find('IdList').findAll('Id')
    ]
    write_ids(ids, output[0])

rule sra_query:
  output:
    "output/{tax_id}/sra/{sra_accession}/query.xml"
  run:
    soup = efetch('sra', wildcards.sra_accession)
    write_soup(soup, output[0])

rule sra_query_json:
  input:
    rules.sra_query.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/query.json"
  run:
    xml2json(input[0], output[0])

rule sra_pluck_biosample:
  input:
    rules.sra_query.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/biosample_id.txt"
  run:
    soup = read_soup(input[0])
    tag = soup.find('EXTERNAL_ID', {'namespace': 'BioSample'})
    id_ = tag.text.strip()
    write_ids([id_], output[0])

rule sra_biosample_fetch:
  input:
    rules.sra_pluck_biosample.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/biosample.xml"
  run:
    with open(input[0]) as biosample_id_file:
      biosample_id = biosample_id_file.read()
    soup = efetch('biosample', biosample_id)
    write_soup(soup, output[0])

rule sra_biosample_json:
  input:
    rules.sra_biosample_fetch.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/biosample.json"
  run:
    xml2json(input[0], output[0])
    

def fetch_instrument_from_sra(sra_query):
  platform = sra_query['EXPERIMENT']['PLATFORM']
  key = list(platform.keys())[0]
  instrument = key + ' - ' + platform[key]['INSTRUMENT_MODEL']
  return instrument


def convert_sample_attributes_to_dict(sra_biosample):
  attribute_list = sra_biosample['Attributes']['Attribute']
  return {
    attribute['@attribute_name']: attribute['#text']
    for attribute in attribute_list
  }


def safe_fetch(dictionary, keys):
  if type(keys) == str:
    return '-' if not keys in dictionary else dictionary[keys]
  else:
    for key in keys:
      if key in dictionary:
        return dictionary[key]
    return '-'

def scrape_sra_query(sra_query):
  sra_query = sra_query['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
  return {
    'sra_run_id': sra_query['RUN_SET']['RUN']['@accession'],
    'instrument': fetch_instrument_from_sra(sra_query),
    'bioproject': sra_query['STUDY']['@alias'],
    'sample_name': sra_query['RUN_SET']['RUN']['Pool']['Member']['@sample_name'],
    'collected_by': sra_query['SUBMISSION']['@center_name']
  }


def scrape_sra_biosample(sra_biosample):
  sra_biosample = sra_biosample['BioSampleSet']['BioSample']
  sample_dict = convert_sample_attributes_to_dict(sra_biosample)
  return {
    'organism_name': sra_biosample['Description']['Organism']['OrganismName'],
    'biosample': sra_biosample['@accession'],
    'strain': safe_fetch(sample_dict, 'strain'),
    'isolate': safe_fetch(sample_dict, 'isolate'),
    'isolation_source': safe_fetch(sample_dict, 'isolation_source'),
    'collection_date': safe_fetch(sample_dict, 'collection_date'),
    'geo_loc_name': safe_fetch(sample_dict, 'geo_loc_name'),
    'isolation_source': safe_fetch(sample_dict, 'isolation_source'),
    'lat_lon': safe_fetch(sample_dict, 'lat_lon'),
    'culture_collection': safe_fetch(sample_dict, 'culture_collection'),
    'host': safe_fetch(sample_dict, ['lab_host', 'host']),
    'host_age': safe_fetch(sample_dict, 'host_age'),
    'host_description': safe_fetch(sample_dict, 'host_description'),
    'host_disease': safe_fetch(sample_dict, ['host_disease', 'disease']),
    'host_disease_outcome': safe_fetch(sample_dict, 'host_disease_outcome'),
    'host_disease_stage': safe_fetch(sample_dict, 'host_disease_stage'),
    'host_health_state': safe_fetch(sample_dict, 'host_health_state'),
    'host_sex': safe_fetch(sample_dict, 'host_sex'),
    'id_method': safe_fetch(sample_dict, 'identification_method')
  }


rule biosample_meta_row:
  input:
    sra=rules.sra_query_json.output[0],
    biosample=rules.sra_biosample_json.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/row.json"
  run:
    sra = read_json(input.sra)
    biosample = read_json(input.biosample)
    scraped = scrape_sra_query(sra)
    scraped.update(scrape_sra_biosample(biosample))
    scraped.update({
      'taxonomy_id': wildcards.tax_id,
      'schema_version': 'v0.8',
      'lineage': lineage[int(wildcards.tax_id)],
      'genome_assembly_id': '-'
    })
    write_json(scraped, output[0])


def build_table(header, rows, out):
    with open(header) as f:
      fieldnames = f.read().strip().split('\t')
    tsv_file = open(out, 'w')
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for row_filepath in rows:
      with open(row_filepath) as json_file:
        row_data = json.load(json_file)
      writer.writerow(row_data)
    tsv_file.close()


rule marburg_sra_biosample_table:
  input:
    header="input/v0.8_biosampleMeta_header.tsv",
    rows=expand(
      "output/11269/sra/{sra_accession}/row.json",
      sra_accession=read_ids('input/sra_accessions/11269.txt')
    )
  output:
    "output/11269/sra/biosampleMeta_PL.tsv"
  run:
    build_table(input.header, input.rows, output[0])

rule influenzaA_sra_biosample_table:
  input:
    header="input/v0.8_biosampleMeta_header.tsv",
    rows=expand(
      "output/11320/sra/{sra_accession}/row.json",
      sra_accession=read_ids('input/sra_accessions/11320.txt')
    )
  output:
    "output/11320/sra/biosampleMeta_PL.tsv"
  run:
    build_table(input.header, input.rows, output[0])

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
