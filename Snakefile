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
  bp_accession="[^/]+",
  tax_id="[^/]+"


lineage = {
  11269: 'Marburg marburgvirus|Marburgvirus|Filoviridae|Mononegavirales|Monjiviricetes|Haploviricotina|Negarnaviricota|Orthornavirae|Riboviria|Viruses',
  11320: 'Influenza A virus|Alphainfluenzavirus|Orthomyxoviridae|Articulavirales|Insthoviricetes|Polyploviricotina|Negarnaviricota|Orthornavirae|Riboviria|Viruses',
  433733: 'human lung metagenome',
  2697049: 'Viruses|Riboviria|Orthornavirae|Pisuviricota|Pisoniviricetes|Nidovirales|Cornidovirineae|Coronaviridae|Orthocoronavirinae|Betacoronavirus|Sarbecovirus|Severe acute respiratory syndrome-related coronavirus',
  211044: 'Influenza A virus|Alphainfluenzavirus|Orthomyxoviridae|Articulavirales|Insthoviricetes|Polyploviricotina|Negarnaviricota|Orthornavirae|Riboviria|Viruses'
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


def esearch(db, term, retmax=1000000):
  handle = Entrez.esearch(db=db, term=term, retmax=retmax)
  raw_xml = handle.read()
  handle.close()
  soup = BeautifulSoup(raw_xml, 'xml')
  return soup


def esummary(db, accession):
  handle = Entrez.esummary(db=db, id=accession, report="full")
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
  if not xml:
    write_ids(['[]'], json_filename)
    return
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


def get_key_from_attribute(attribute):
  if '@harmonized_name' in attribute:
    return attribute['@harmonized_name']
  return attribute['@attribute_name']


def convert_sample_attributes_to_dict(sra_biosample):
  if sra_biosample['Attributes'] is None:
    return {}
  attribute_list = sra_biosample['Attributes']['Attribute']
  if type(attribute_list) == dict:
    attribute_list = [attribute_list]
  return {
    get_key_from_attribute(attribute): attribute['#text']
    for attribute in attribute_list
  }


def scrape_biosample(sra_biosample):
  if sra_biosample == []:
    organism_name = '-'
    biosample = '-'
    sample_dict = {}
  else:
    if sra_biosample['BioSampleSet'] is None:
      organism_name = '-'
      biosample = '-'
      sample_dict = {}
    else:
      base = sra_biosample['BioSampleSet']['BioSample']
      sample_dict = convert_sample_attributes_to_dict(base)
      organism_name = safe_fetch(
        base['Description']['Organism'],
        ['OrganismName', '@taxonomy_name']
      )
      biosample = base['@accession']

  return {
    'organism_name': organism_name,
    'biosample': biosample,
    'strain': safe_fetch(sample_dict, 'strain'),
    'isolate': safe_fetch(sample_dict, 'isolate'),
    'isolation_source': safe_fetch(sample_dict, 'isolation_source'),
    'collection_date': safe_fetch(sample_dict, ['collection_date', 'collection date']),
    'geo_loc_name': safe_fetch(sample_dict, ['geo_loc_name', 'geographic location (country and/or sea)']),
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
    'id_method': safe_fetch(sample_dict, 'identification_method'),
    'sample_name': safe_fetch(sample_dict, 'sample_name')
  }


def scrape_assembly(assembly):
  base = assembly['eSummaryResult']['DocumentSummarySet']['DocumentSummary']
  gb_key = 'GB_BioProjects'
  rs_key = 'RS_BioProjects'
  has_gbbp_key = gb_key in base and base[gb_key] != None
  has_rsbp_key = rs_key in base and base[rs_key] != None
  if has_gbbp_key or has_rsbp_key:
    key = gb_key if has_gbbp_key else rs_key
    bioproject = base[key]['Bioproj']['BioprojectAccn']
  else:
    bioproject = '-'
  has_organism = 'Organism' in base
  organism = '-' if not has_organism else base['Organism']
  return {
    'bioproject': bioproject,
    'organism_name': safe_fetch(base, 'Organism'),
    'genome_assembly_id': safe_fetch(base, 'AssemblyAccession'),
    'instrument': '-',
    'collected_by': '-'
  }


rule db_search:
  output:
    "db/{database}/{id_}/search.xml"
  run:
    soup = esearch(wildcards.database, wildcards.id_)
    write_soup(soup, output[0])

rule db_search_json:
  input:
    rules.db_search.output[0]
  output:
    "db/{database}/{id_}/search.json"
  run:
    xml2json(input[0], output[0]) 

rule db_fetch:
  output:
    "db/{database}/{id_}/fetch.xml"
  run:
    soup = efetch(wildcards.database, wildcards.id_)
    write_soup(soup, output[0])

rule db_fetch_json:
  input:
    rules.db_fetch.output[0]
  output:
    "db/{database}/{id_}/fetch.json"
  run:
    xml2json(input[0], output[0]) 

rule taxon_search:
  output:
    "output/{tax_id}/{database}/search.xml"
  run:
    tax_id = 'txid' + wildcards.tax_id + '[ORGN]'
    soup = esearch(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_search_json:
  input:
    rules.taxon_search.output[0]
  output:
    "output/{tax_id}/{database}/search.json"
  run:
    xml2json(input[0], output[0])

rule taxon_summary:
  output:
    "output/{tax_id}/{database}/summary.xml"
  run:
    tax_id = wildcards.tax_id
    soup = esummary(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_summary_json:
  input:
    rules.taxon_summary.output[0]
  output:
    "output/{tax_id}/{database}/summary.json"
  run:
    xml2json(input[0], output[0])

rule taxon_fetch:
  output:
    "output/{tax_id}/{database}/fetch.xml"
  run:
    tax_id = wildcards.tax_id
    soup = efetch(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_fetch_json:
  input:
    rules.taxon_fetch.output[0]
  output:
    "output/{tax_id}/{database}/fetch.json"
  run:
    xml2json(input[0], output[0])

checkpoint accessions:
  input:
    rules.taxon_search_json.output[0]
  output:
    "output/{tax_id}/{database}/accessions.txt"
  run:
    query = read_json(input[0])
    id_list = query['eSearchResult']['IdList']['Id']
    write_ids(id_list, output[0])

rule assembly_query:
  output:
    "output/{tax_id}/assembly/{assembly_id}/query.xml"
  run:
    soup = esummary('assembly', wildcards.assembly_id)
    write_soup(soup, output[0])

rule assembly_query_json:
  input:
    rules.assembly_query.output[0]
  output:
    "output/{tax_id}/assembly/{assembly_id}/query.json"
  run:
    xml2json(input[0], output[0])

rule assembly_pluck_biosample_id:
  input:
    rules.assembly_query_json.output[0]
  output:
    "output/{tax_id}/assembly/{assembly_id}/biosample_id.txt"
  run:
    query = read_json(input[0])
    base = query['eSummaryResult']['DocumentSummarySet']['DocumentSummary']
    biosample_id = '' if not 'BioSampleId' in base else base['BioSampleId']
    if biosample_id is None:
      biosample_id = ''
    write_ids([biosample_id], output[0])

rule assembly_biosample_fetch:
  input:
    rules.assembly_pluck_biosample_id.output[0]
  output:
    "output/{tax_id}/assembly/{assembly_id}/biosample.xml"
  run:
    with open(input[0]) as biosample_id_file:
      biosample_id = biosample_id_file.read()
    if biosample_id:
      soup = efetch('biosample', biosample_id)
      write_soup(soup, output[0])
    else:
      write_ids([biosample_id], output[0])

rule assembly_biosample_json:
  input:
    rules.assembly_biosample_fetch.output[0]
  output:
    "output/{tax_id}/assembly/{assembly_id}/biosample.json"
  run:
    xml2json(input[0], output[0])

rule assembly_row:
  input:
    assembly=rules.assembly_query_json.output[0],
    biosample=rules.assembly_biosample_json.output[0]
  output:
    "output/{tax_id}/assembly/{assembly_id}/row.json"
  run:
    assembly = read_json(input.assembly)
    biosample = read_json(input.biosample)
    scraped = scrape_biosample(biosample)
    scraped.update(scrape_assembly(assembly))
    scraped.update({
      'taxonomy_id': wildcards.tax_id,
      'schema_version': 'v0.8',
      'lineage': lineage[int(wildcards.tax_id)],
      'sra_run_id': '-',
    })
    write_json(scraped, output[0])

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


def safe_fetch(dictionary, keys):
  if type(keys) == str:
    return '-' if not keys in dictionary else dictionary[keys]
  else:
    for key in keys:
      if key in dictionary:
        return dictionary[key]
    return '-'


def deep_safe_fetch(dictionary, keys):
  for key in keys:
    if key in dictionary:
      dictionary = dictionary[key]
      if type(dictionary) != dict:
        return dictionary
    else:
      return '-'

def scrape_sra_query(sra_query, sra_accession):
  sra_query = sra_query['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
  run_ = sra_query['RUN_SET']['RUN']
  if type(run_) == list:
    run_ = [r for r in run_ if r['@accession'] == sra_accession][0]
  return {
    'instrument': fetch_instrument_from_sra(sra_query),
    'bioproject': sra_query['STUDY']['@alias'],
    'sample_name':  deep_safe_fetch(run_, ['Pool', 'Member' ,'@sample_name']),
    'collected_by': deep_safe_fetch(sra_query, ['SUBMISSION', '@center_name'])
  }

rule biosample_meta_row:
  input:
    sra=rules.sra_query_json.output[0],
    biosample=rules.sra_biosample_json.output[0]
  output:
    "output/{tax_id}/sra/{sra_accession}/row.json"
  run:
    sra_query = read_json(input.sra)
    biosample = read_json(input.biosample)
    scraped = scrape_sra_query(sra_query, wildcards.sra_accession)
    scraped.update(scrape_biosample(biosample))
    scraped.update({
      'sra_run_id': wildcards.sra_accession,
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


def table_input(wildcards):
  parameters = (wildcards.tax_id, wildcards.database)
  input_accessions = checkpoints.accessions.get(**wildcards).output[0]
  files = expand(
    "output/%s/%s/{accession}/row.json" % parameters,
    accession=read_ids(input_accessions)
  )
  return files


rule sub_table:
  input:
    header="input/v0.8_biosampleMeta_header.tsv",
    accessions="output/{tax_id}/{database}/accessions.txt",
    rows=table_input
  output:
    "output/{tax_id}/{database}/biosampleMeta_PL.tsv"
  run:
    build_table(input.header, input.rows, output[0])

#rule marburg_assembly_table:
#  input:
#    header="input/v0.8_biosampleMeta_header.tsv",
#    rows=expand(
#      "output/11269/assembly/{assembly_accession}/row.json",
#      assembly_accession=read_ids('output/11269/assembly/accessions.txt')
#    )
#  output:
#    "output/11269/assembly/biosampleMeta_PL.tsv"
#  run:
#    build_table(input.header, input.rows, output[0])
#
#rule influenzaA_sra_table:
#  input:
#    header="input/v0.8_biosampleMeta_header.tsv",
#    rows=expand(
#      "output/11320/sra/{sra_accession}/row.json",
#      sra_accession=read_ids('input/sra_accessions/11320.txt')
#    )
#  output:
#    "output/11320/sra/biosampleMeta_PL.tsv"
#  run:
#    build_table(input.header, input.rows, output[0])
#
#rule influenzaA_assembly_table:
#  input:
#    header="input/v0.8_biosampleMeta_header.tsv",
#    rows=expand(
#      "output/11320/assembly/{assembly_accession}/row.json",
#      assembly_accession=read_ids('output/11320/assembly/accessions.txt')
#    )
#  output:
#    "output/11320/assembly/biosampleMeta_PL.tsv"
#  run:
#    build_table(input.header, input.rows, output[0])
#
#rule master_table:
#  input:
#    header="input/v0.8_biosampleMeta_header.tsv",
#    marburg_assembly=rules.marburg_assembly_table.output[0],
#    influenzaA_assembly=rules.influenzaA_assembly_table.output[0],
#    marburg_sra=rules.marburg_sra_table.output[0],
#    influenzaA_sra=rules.influenzaA_sra_table.output[0]
#  output:
#    "output/biosampleMeta_PL.tsv"
#  shell:
#    """
#      cp {input.header} {output}
#      tail -n +2 {input.marburg_assembly} >> {output}
#      tail -n +2 {input.influenzaA_assembly} >> {output}
#      tail -n +2 {input.marburg_sra} >> {output}
#      tail -n +2 {input.influenzaA_sra} >> {output}
#    """
