import json

from Bio import Entrez
from bs4 import BeautifulSoup
import xmltodict


def efetch(db, accession):
    handle = Entrez.efetch(db=db, id=accession, rettype="xml")
    raw_xml = handle.read()
    handle.close()
    soup = BeautifulSoup(raw_xml, "xml")
    return soup


def elink(db, dbfrom, id_, retmax=1000):
    handle = Entrez.elink(db=db, dbfrom=dbfrom, id=id_, retmax=retmax)
    raw_xml = handle.read()
    handle.close()
    soup = BeautifulSoup(raw_xml, "xml")
    return soup


def esearch(db, term, retmax=1000000):
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    raw_xml = handle.read()
    handle.close()
    soup = BeautifulSoup(raw_xml, "xml")
    return soup


def esummary(db, accession):
    handle = Entrez.esummary(db=db, id=accession, report="full")
    raw_xml = handle.read()
    handle.close()
    soup = BeautifulSoup(raw_xml, "xml")
    return soup


def read_soup(filename):
    with open(filename) as f:
        soup = BeautifulSoup(f, "xml")
    return soup


def write_soup(soup, filename):
    pretty = soup.prettify()
    with open(filename, "w") as f:
        f.write(pretty)


def read_ids(filename):
    with open(filename) as f:
        ids = f.read().splitlines()
    return ids


def write_ids(ids, filename):
    with open(filename, "w") as f:
        f.write("\n".join(ids))


def read_json(filename):
    with open(filename) as json_file:
        result = json.load(json_file)
    return result


def write_json(result, filename):
    with open(filename, "w") as json_file:
        json.dump(result, json_file, indent=2)


def xml2json(xml_filename, json_filename):
    with open(xml_filename) as xml_file:
        xml = xml_file.read()
    if not xml:
        write_ids(["[]"], json_filename)
        return
    converted = xmltodict.parse(xml)
    with open(json_filename, "w") as json_file:
        json.dump(converted, json_file, indent=2)


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


def fetch_instrument_from_sra(sra_query):
  platform = sra_query['EXPERIMENT']['PLATFORM']
  key = list(platform.keys())[0]
  instrument = key + ' - ' + platform[key]['INSTRUMENT_MODEL']
  return instrument


def scrape_sra_query(sra_query, sra_accession=None):
  sra_query = sra_query['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
  run_ = sra_query['RUN_SET']['RUN']
  if type(run_) == list:
    run_ = [r for r in run_ if r['@accession'] == sra_accession][0]
  return {
    'instrument': fetch_instrument_from_sra(sra_query),
    'bioproject': sra_query['STUDY']['IDENTIFIERS']['EXTERNAL_ID']['#text'],
    'biosample': pluck_biosample_from_sra(sra_query),
    'sample_name':  deep_safe_fetch(run_, ['Pool', 'Member' ,'@sample_name']),
    'organism': sra_query['Pool']['Member']['@organism'],
    'tax_id': sra_query['SAMPLE']['SAMPLE_NAME']['TAXON_ID'],
    'collected_by': deep_safe_fetch(sra_query, ['SUBMISSION', '@center_name'])
  }


def pluck_biosample_from_sra(input_sra):
    id_ = input_sra["SAMPLE"]["IDENTIFIERS"]["EXTERNAL_ID"]["#text"]
    return id_


def get_key_from_attribute(attribute):
    if "@harmonized_name" in attribute:
        return attribute["@harmonized_name"]
    return attribute["@attribute_name"]


def convert_sample_attributes_to_dict(sra_biosample):
    if sra_biosample["Attributes"] is None:
        return {}
    attribute_list = sra_biosample["Attributes"]["Attribute"]
    if type(attribute_list) == dict:
        attribute_list = [attribute_list]
    return {
        get_key_from_attribute(attribute): attribute["#text"]
        for attribute in attribute_list
    }


def scrape_biosample(sra_biosample):
    if sra_biosample == []:
        organism_name = "-"
        biosample = "-"
        sample_dict = {}
    else:
        if sra_biosample["BioSampleSet"] is None:
            organism_name = "-"
            biosample = "-"
            sample_dict = {}
        else:
            base = sra_biosample["BioSampleSet"]["BioSample"]
            sample_dict = convert_sample_attributes_to_dict(base)
            organism_name = safe_fetch(
                base["Description"]["Organism"], ["OrganismName", "@taxonomy_name"]
            )
            biosample = base["@accession"]

    return {
        "organism_name": organism_name,
        "biosample": biosample,
        "strain": safe_fetch(sample_dict, "strain"),
        "isolate": safe_fetch(sample_dict, "isolate"),
        "isolation_source": safe_fetch(sample_dict, "isolation_source"),
        "collection_date": safe_fetch(
            sample_dict, ["collection_date", "collection date"]
        ),
        "geo_loc_name": safe_fetch(
            sample_dict, ["geo_loc_name", "geographic location (country and/or sea)"]
        ),
        "isolation_source": safe_fetch(sample_dict, "isolation_source"),
        "lat_lon": safe_fetch(sample_dict, "lat_lon"),
        "culture_collection": safe_fetch(sample_dict, "culture_collection"),
        "host": safe_fetch(sample_dict, ["lab_host", "host"]),
        "host_age": safe_fetch(sample_dict, "host_age"),
        "host_description": safe_fetch(sample_dict, "host_description"),
        "host_disease": safe_fetch(sample_dict, ["host_disease", "disease"]),
        "host_disease_outcome": safe_fetch(sample_dict, "host_disease_outcome"),
        "host_disease_stage": safe_fetch(sample_dict, "host_disease_stage"),
        "host_health_state": safe_fetch(sample_dict, "host_health_state"),
        "host_sex": safe_fetch(sample_dict, "host_sex"),
        "id_method": safe_fetch(sample_dict, "identification_method"),
        "sample_name": safe_fetch(sample_dict, "sample_name"),
    }


def scrape_assembly(assembly):
    base = assembly["eSummaryResult"]["DocumentSummarySet"]["DocumentSummary"]
    gb_key = "GB_BioProjects"
    rs_key = "RS_BioProjects"
    has_gbbp_key = gb_key in base and base[gb_key] != None
    has_rsbp_key = rs_key in base and base[rs_key] != None
    if has_gbbp_key or has_rsbp_key:
        key = gb_key if has_gbbp_key else rs_key
        bioproject = base[key]["Bioproj"]["BioprojectAccn"]
    else:
        bioproject = "-"
    has_organism = "Organism" in base
    organism = "-" if not has_organism else base["Organism"]
    return {
        "bioproject": bioproject,
        "organism_name": safe_fetch(base, "Organism"),
        "genome_assembly_id": safe_fetch(base, "AssemblyAccession"),
        "instrument": "-",
        "collected_by": "-",
    }

def extract_xml_link_ids(soup):
    return [
      id_.text.strip()
      for id_ in soup.find('LinkSetDb').findAll('Id')
    ]


def extract_json_link_ids(link_dict):
    link_set = link_dict['eLinkResult']['LinkSet']
    if not 'LinkSetDb' in link_set:
        return []
    links = link_set['LinkSetDb']['Link']
    if type(links) == list:
        link_list = [v['Id'] for v in links]
        return link_list
    return [links['Id']]
