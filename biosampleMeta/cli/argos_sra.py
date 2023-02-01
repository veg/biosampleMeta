

def argos_sra_cli():
    parser = argparse.ArgumentParser(description='Scrape SRA data for ArgosDB.')
    parser.add_argument('--input', type=str, help='input JSON')
    parser.add_argument('--output', type=str, help='output JSON')
    args = parser.parse_args()
    with open(args.input) as json_file:
        argos = scrape_sra(json.load(json_file))
    with open(args.output, 'w') as json_file:
        json.dump(argos, json_file, indent=2)



