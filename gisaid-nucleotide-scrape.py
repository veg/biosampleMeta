import json
import xmltodict


nucleotide_file = open('input/gisaid-nucleotide.xml')

lines_to_keep = []
sets_kept = 0
tag_counter = 0
tag = 'Seq-entry_set'
in_tag = False
process_chunk = False
for line in nucleotide_file:
    if '<%s>' % tag in line:
        in_tag = True
        tag_counter += 1
    if in_tag:
        lines_to_keep.append(line)
    if '</%s>' % tag in line:
        in_tag = False
        tag_counter -= 1
        if tag_counter == 0:
            process_chunk = True
    if process_chunk:
        process_chunk = False
        xml = ''.join(lines_to_keep)
        xml_filename = 'output/nucleotide-parsing/set-%d.xml' % sets_kept
        with open(xml_filename, 'w') as xml_file:
            xml_file.write(xml)
        json_filename = 'output/nucleotide-parsing/set-%d.json' % sets_kept
        with open(json_filename, 'w') as json_file:
            as_dict = xmltodict.parse(xml)
            json_file.write(json.dumps(as_dict))
        sets_kept += 1
        print('Wrote %d lines to %s...' % (len(lines_to_keep), json_filename))
        tag_counter = 0
        lines_to_keep = []
nucleotide_file.close()
