import json

from Bio import Entrez
from bs4 import BeautifulSoup
import xmltodict


def efetch(db, accession):
    handle = Entrez.efetch(db=db, id=accession)
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
