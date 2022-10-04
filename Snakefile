import os
import csv

from Bio import Entrez

from biosampleMeta import *


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

rule bioproject_query:
  output:
    "output/tax_id/{tax_id}/bioprojects/search.xml",
  run:
    tax_id = 'txid' + wildcards.tax_id
    soup = esearch('bioproject', tax_id)
    write_soup(soup, output[0])

rule bioproject_ids_from_query:
  input:
    rules.bioproject_query.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/ids.txt"
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
      "output/tax_id/{tax_id}/bioprojects/ids.txt",
      tax_id=read_ids('input/tax_ids.txt')
    )

rule bioproject_xml_fetch:
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/bioproject.xml"
  run:
    soup = efetch('bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_query:
  input:
    rules.bioproject_xml_fetch.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/sra_accessions.xml"
  run:
    soup = elink('sra', 'bioproject', id_=wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_sra_ids:
  input:
    rules.bioproject_sra_query.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/sra_accessions.txt"
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
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/biosample_links.xml"
  run:
    soup = elink('biosample', 'bioproject', wildcards.bp_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_ids:
  input:
    rules.bioproject_biosample_query.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/biosample_accessions.txt"
  run:
    soup = read_soup(input[0])
    samn_ids = [
      id_.text.strip()
      for id_ in soup.find('LinkSetDb').findAll('Id')
    ]
    write_ids(samn_ids, output[0])
 
rule bioproject_biosample_entry:
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/entry.xml"
  run:
    soup = efetch('biosample', wildcards.bs_accession)
    write_soup(soup, output[0])

rule bioproject_biosample_geo_accession:
  input:
    rules.bioproject_biosample_entry.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/geo_accession.txt"
  run:
    soup = read_soup(input[0])
    geo_accession = soup.find('Id', {'db': "GEO"}).text.strip()
    write_ids([geo_accession], output[0])

rule bioproject_biosample_geo_xml:
  input:
    rules.bioproject_biosample_geo_accession.output[0]
  output:
    "output/tax_id/{tax_id}/bioprojects/{bp_accession}/biosamples/{bs_accession}/geo.xml"
  shell:
    """
      sleep 3;
      GEO_ACCESSION=$(cat {input})
      wget "https://www.ncbi.nlm.nih.gov/geo/tools/geometa.cgi?acc=$GEO_ACCESSION&scope=full&mode=miniml" -O {output}
    """

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
    "output/tax_id/{tax_id}/{database}/search.xml"
  run:
    tax_id = 'txid' + wildcards.tax_id + '[ORGN]'
    soup = esearch(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_search_json:
  input:
    rules.taxon_search.output[0]
  output:
    "output/tax_id/{tax_id}/{database}/search.json"
  run:
    xml2json(input[0], output[0])

rule taxon_summary:
  output:
    "output/tax_id/{tax_id}/{database}/summary.xml"
  run:
    tax_id = wildcards.tax_id
    soup = esummary(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_summary_json:
  input:
    rules.taxon_summary.output[0]
  output:
    "output/tax_id/{tax_id}/{database}/summary.json"
  run:
    xml2json(input[0], output[0])

rule taxon_fetch:
  output:
    "output/tax_id/{tax_id}/{database}/fetch.xml"
  run:
    tax_id = wildcards.tax_id
    soup = efetch(wildcards.database, tax_id)
    write_soup(soup, output[0])

rule taxon_fetch_json:
  input:
    rules.taxon_fetch.output[0]
  output:
    "output/tax_id/{tax_id}/{database}/fetch.json"
  run:
    xml2json(input[0], output[0])

checkpoint accessions:
  input:
    rules.taxon_search_json.output[0]
  output:
    "output/tax_id/{tax_id}/{database}/accessions.txt"
  run:
    query = read_json(input[0])
    id_list = query['eSearchResult']['IdList']['Id']
    write_ids(id_list, output[0])

rule assembly_query:
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/query.xml"
  run:
    soup = esummary('assembly', wildcards.assembly_id)
    write_soup(soup, output[0])

rule assembly_query_json:
  input:
    rules.assembly_query.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/query.json"
  run:
    xml2json(input[0], output[0])

rule assembly_pluck_biosample_id:
  input:
    rules.assembly_query_json.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/biosample_id.txt"
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
    "output/tax_id/{tax_id}/assembly/{assembly_id}/biosample.xml"
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
    "output/tax_id/{tax_id}/assembly/{assembly_id}/biosample.json"
  run:
    xml2json(input[0], output[0])

rule assembly_pluck_sraid_from_biosample:
  input:
    rules.assembly_biosample_json.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra_id.txt"
  run:
    biosample = read_json(input[0])
    keys = ['BioSampleSet', 'BioSample', 'Ids', 'Id']
    ids = deep_safe_fetch(biosample, keys)
    id_ = [an_id['#text'] for an_id in ids if an_id['@db'] == 'SRA']
    write_ids(id_, output[0])

rule assembly_sra_fetch:
  input:
    rules.assembly_pluck_sraid_from_biosample.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra.xml"
  run:
    id_ = read_ids(input[0])[0]
    print(id_)
    soup = efetch('sra', id_)
    write_soup(soup, output[0])

rule assembly_sra_query:
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra_query.xml"
  run:
    soup = elink('sra', 'assembly', id_=wildcards.assembly_id)
    write_soup(soup, output[0])

rule assembly_sra_query_json:
  input:
    rules.assembly_sra_query.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra_query.json"
  run:
    xml2json(input[0], output[0])

rule assembly_row:
  input:
    assembly=rules.assembly_query_json.output[0],
    biosample=rules.assembly_biosample_json.output[0],
    sra=rules.assembly_sra_query_json.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/row.json"
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

rule tax_id_sra_query:
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/query.xml"
  run:
    soup = efetch('sra', wildcards.sra_accession)
    write_soup(soup, output[0])

rule tax_id_sra_query_json:
  input:
    rules.tax_id_sra_query.output[0]
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/query.json"
  run:
    xml2json(input[0], output[0])

rule sra_pluck_biosample:
  input:
    rules.tax_id_sra_query.output[0]
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/biosample_id.txt"
  run:
    soup = read_soup(input[0])
    tag = soup.find('EXTERNAL_ID', {'namespace': 'BioSample'})
    id_ = tag.text.strip()
    write_ids([id_], output[0])

rule sra_biosample_fetch:
  input:
    rules.sra_pluck_biosample.output[0]
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/biosample.xml"
  run:
    with open(input[0]) as biosample_id_file:
      biosample_id = biosample_id_file.read()
    soup = efetch('biosample', biosample_id)
    write_soup(soup, output[0])

rule sra_biosample_json:
  input:
    rules.sra_biosample_fetch.output[0]
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/biosample.json"
  run:
    xml2json(input[0], output[0])

rule biosample_meta_row:
  input:
    sra=rules.tax_id_sra_query_json.output[0],
    biosample=rules.sra_biosample_json.output[0]
  output:
    "output/tax_id/{tax_id}/sra/{sra_accession}/row.json"
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
    "output/tax_id/%s/%s/{accession}/row.json" % parameters,
    accession=read_ids(input_accessions)
  )
  return files


rule sub_table:
  input:
    header="input/v0.8_biosampleMeta_header.tsv",
    accessions="output/tax_id/{tax_id}/{database}/accessions.txt",
    rows=table_input
  output:
    "output/tax_id/{tax_id}/{database}/biosampleMeta_PL.tsv"
  run:
    build_table(input.header, input.rows, output[0])

rule sra_query:
  output:
    "output/sra/{sra_accession}/query.xml"
  run:
    soup = efetch('sra', wildcards.sra_accession)
    write_soup(soup, output[0])

rule sra_query_json:
  input:
    rules.sra_query.output[0]
  output:
    "output/sra/{sra_accession}/query.json"
  run:
    xml2json(input[0], output[0])

rule sra_pluck_biosample_id:
  input:
    rules.sra_query_json.output[0]
  output:
    "output/sra/{sra_accession}/biosample_id.txt"
  run:
    pluck_biosample_from_sra(input[0], output[0])

rule sra_biosample_table:
  input:
    expand(
      "output/sra/{sra_accession}/biosample_id.txt",
      sra_accession=read_ids('input/accessions.txt')
    )
