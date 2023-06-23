import os
import csv

from Bio import Entrez

from biosampleMeta import *


Entrez.email = os.environ.get('ENTREZ_EMAIL') or None
Entrez.api_key = os.environ.get('ENTREZ_API_KEY') or None


wildcard_constraints:
  bp_accession="[^/]+",
  tax_id="[^/]+",
  bioproject_id="[^/]+",
  database="[^/]+",
  sra_accession="[^/]+",
  db_id="[^/]+",
  otherdb="[^/]+",
  otherdb_id="[^/]+"

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
    ids = extract_xml_link_ids(soup)
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
    "output/db/{database}/{id_}/search.xml"
  run:
    soup = esearch(wildcards.database, wildcards.id_)
    write_soup(soup, output[0])

rule db_search_json:
  input:
    rules.db_search.output[0]
  output:
    "output/db/{database}/{id_}/search.json"
  run:
    xml2json(input[0], output[0]) 

rule db_fetch:
  output:
    "output/db/{database}/{id_}/fetch.xml"
  run:
    soup = efetch(wildcards.database, wildcards.id_)
    write_soup(soup, output[0])

rule db_fetch_json:
  input:
    rules.db_fetch.output[0]
  output:
    "output/db/{database}/{id_}/fetch.json"
  run:
    xml2json(input[0], output[0]) 

rule db_link:
  output:
    xml="output/db/{database}/{id_}/{database_from}/link.xml",
    json="output/db/{database}/{id_}/{database_from}/link.json"
  run:
    soup = elink(wildcards.database_from, wildcards.database, wildcards.id_)
    write_soup(soup, output.xml)
    xml2json(output.xml, output.json)

rule db_summary:
  output:
    "output/db/{database}/{id_}/summary.xml"
  run:
    soup = esummary(wildcards.database, wildcards.id_)
    write_soup(soup, output[0])

rule db_summary_json:
  input:
    rules.db_summary.output[0]
  output:
    "output/db/{database}/{id_}/summary.json"
  run:
    xml2json(input[0], output[0])
#
#checkpoint accessions:
#  input:
#    rules.taxon_search_json.output[0]
#  output:
#    "output/tax_id/{tax_id}/{database}/accessions.txt"
#  run:
#    query = read_json(input[0])
#    id_list = query['eSearchResult']['IdList']['Id']
#    write_ids(id_list, output[0])
#

rule assembly_pluck_biosample_id:
  input:
    "output/db/assembly/{assembly_id}/summary.json"
  output:
    "output/db/assembly/{assembly_id}/biosample_id.txt"
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
    id_ = [
      an_id['#text']
      for an_id in ids
      if safe_fetch(an_id, '@db') == 'SRA'
    ]
    write_ids(id_, output[0])

rule assembly_sra_fetch:
  input:
    rules.assembly_pluck_sraid_from_biosample.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra.xml"
  run:
    id_ = read_ids(input[0])[0].strip()
    soup = esearch('sra', id_)
    write_soup(soup, output[0])

rule assembly_sra_fetch_json:
  input:
    rules.assembly_sra_fetch.output[0]
  output:
    "output/tax_id/{tax_id}/assembly/{assembly_id}/sra.json"
  run:
    xml2json(input[0], output[0])

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

#rule assembly_row:
#  input:
#    assembly=rules.assembly_query_json.output[0],
#    biosample=rules.assembly_biosample_json.output[0],
#    sra=rules.assembly_sra_query_json.output[0]
#  output:
#    "output/{db}/assembly/{assembly_id}/row.json"
#  run:
#    assembly = read_json(input.assembly)
#    biosample = read_json(input.biosample)
#    scraped = scrape_biosample(biosample)
#    scraped.update(scrape_assembly(assembly))
#    scraped.update({
#      'taxonomy_id': wildcards.tax_id,
#      'schema_version': 'v0.8',
#      'lineage': lineage[int(wildcards.tax_id)],
#      'sra_run_id': '-',
#    })
#    write_json(scraped, output[0])

rule bioproject_db_links_xml:
  output:
    "output/db/bioproject/{bioproject_id}/{db}/links.xml"
  run:
    soup = elink(wildcards.db, 'bioproject', id_=wildcards.bioproject_id)
    write_soup(soup, output[0])

rule bioproject_db_links_text:
  input:
    rules.bioproject_db_links_xml.output[0]
  output:
    "output/bioproject/{bioproject_id}/{db}/links.txt"
  run:
    soup = read_soup(input[0])
    ids = extract_xml_link_ids(soup)
    write_ids(ids, output[0])

rule bioproject_db_otherdb_links_xml:
  output:
    "output/bioproject/{bioproject_id}/{db}/{db_id}/{otherdb}/links.xml"
  run:
    soup = elink(wildcards.otherdb, wildcards.db, id_=wildcards.db_id)
    write_soup(soup, output[0])

rule bioproject_db_otherdb_links_text:
  input:
    rules.bioproject_db_otherdb_links_xml.output[0]
  output:
    "output/bioproject/{bioproject_id}/{db}/{db_id}/{otherdb}/links.txt"
  run:
    soup = read_soup(input[0])
    ids = extract_xml_link_ids(soup)
    write_ids(ids, output[0])

#rule bioproject_db_otherdb_link_query:
#  output:
#    "output/bioproject/{bioproject_id}/{db}/{db_id}/{otherdb}/{otherdb_id}/query.xml"
#  run:
#    print(wildcards.otherdb, wildcards.otherdb_id)
#    soup = efetch(wildcards.otherdb, wildcards.otherdb_id)
#    write_soup(soup, output[0])
#
#rule biosample_meta_row:
#  input:
#    sra=rules.tax_id_sra_query_json.output[0],
#    biosample=rules.sra_biosample_json.output[0]
#  output:
#    "output/tax_id/{tax_id}/sra/{sra_accession}/row.json"
#  run:
#    sra_query = read_json(input.sra)
#    biosample = read_json(input.biosample)
#    scraped = scrape_sra_query(sra_query, wildcards.sra_accession)
#    scraped.update(scrape_biosample(biosample))
#    scraped.update({
#      'sra_run_id': wildcards.sra_accession,
#      'taxonomy_id': wildcards.tax_id,
#      'schema_version': 'v0.8',
#      'lineage': lineage[int(wildcards.tax_id)],
#      'genome_assembly_id': '-'
#    })
#    write_json(scraped, output[0])
#
#
#def build_table(header, rows, out):
#    with open(header) as f:
#      fieldnames = f.read().strip().split('\t')
#    y
#    writer.writeheader()
#    for row_filepath in rows:
#      with open(row_filepath) as json_file:
#        row_data = json.load(json_file)
#      writer.writerow(row_data)
#    tsv_file.close()
#
#
#def table_input(wildcards):
#  parameters = (wildcards.tax_id, wildcards.database)
#  input_accessions = checkpoints.accessions.get(**wildcards).output[0]
#  files = expand(
#    "output/tax_id/%s/%s/{accession}/row.json" % parameters,
#    accession=read_ids(input_accessions)
#  )
#  return files
#
#
#rule sub_table:
#  input:
#    header="input/v0.8_biosampleMeta_header.tsv",
#    accessions="output/tax_id/{tax_id}/{database}/accessions.txt",
#    rows=table_input
#  output:
#    "output/tax_id/{tax_id}/{database}/biosampleMeta_PL.tsv"
#  run:
#    build_table(input.header, input.rows, output[0])

rule sra_argos_scrape:
  input:
    rules.db_fetch_json.output[0]
  output:
    "output/db/{database}/{id_}/sra_scrape.json"
  run:
    sra = read_json(input[0])
    scrape = scrape_sra_query(sra)
    write_json(scrape, output[0])

rule sra_pluck_biosample_id:
  input:
    rules.db_fetch_json.output[0]
  output:
    "output/db/{database}/{id_}/biosample_id_from_sra.txt"
  run:
    pluck_biosample_from_sra(input[0], output[0])

rule sra_biosample_table:
  input:
    expand(
      "output/sra/{sra_accession}/biosample_id.txt",
      sra_accession=read_ids('input/accessions.txt')
    )

def argos_biosample_input(wildcards):
  biosample_ids = read_ids("./input/PRJNA231221_biosample_ids.txt")
  nucleotide = [
    "output/db/biosample/%s/nucleotide/link.json" % biosample_id
    for biosample_id in biosample_ids
  ]
  sra = [
    "output/db/biosample/%s/sra/link.json" % biosample_id
    for biosample_id in biosample_ids
  ]
  return nucleotide+sra

rule argos_biosample_links:
  input:
    argos_biosample_input
  output:
    "output/argos_table.tsv"
  run:
    harmonized = {}
    for filename in input:
      id_ = filename.split('/')[3]
      db = filename.split('/')[4]
      link_query = read_json(filename)
      link_list = extract_json_link_ids(link_query)
      link_ids = ', '.join(link_list)
      if db == 'nucleotide':
        harmonized[id_] = {'nucleotide': link_ids}
      else:
        harmonized[id_]['sra'] = link_ids
    tsv_file = open(output[0], 'w')
    fieldnames = ['biosample', 'nucleotide', 'sra']
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for key, value in harmonized.items():
      row = value
      row['biosample'] = key
      writer.writerow(row)
    tsv_file.close()

rule biosample_scrape:
  input:
    "output/db/biosample/{bs_id}/fetch.json"
  output:
    "output/db/biosample/{bs_id}/argos.json"
  run:
    bs = read_json(input[0])
    argos = scrape_biosample(bs)
    write_json(argos, output[0])

rule augment_argos_table:
  input:
    rules.argos_biosample_links.output[0]
  output:
    "output/full_argos_table.tsv"
  run:
    input_tsv_file = open(input[0])
    reader = csv.DictReader(input_tsv_file, delimiter='\t')
    output_tsv_file = open(output[0], 'w')
    fieldnames = ['biosample', 'nucleotide', 'sra', 'organism_name']
    writer = csv.DictWriter(output_tsv_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for row in reader:
      bs_id = row['biosample']
      argos_bs = read_json('output/db/biosample/%s/argos.json' % row['biosample'])
      row.update({'organism_name': argos_bs['organism_name']})
      writer.writerow(row)
