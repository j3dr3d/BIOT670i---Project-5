from Bio import Entrez
import os
import sys
import time
import csv

# Replace with your actual email
Entrez.email = "your.email@example.com"

def get_sra_ids_from_samn(samn_ids, batch_size=200):
    """
    Finds BioSample IDs, links to SRA, and fetches SRR accessions using the
    robust efetch 'runinfo' method.
    """
    srr_accessions = []
    total_ids = len(samn_ids)
    print(f"Starting search for {total_ids} BioSample accessions.")
    
    for i in range(0, total_ids, batch_size):
        batch_samn_ids = samn_ids[i:i + batch_size]
        id_string = " OR ".join(batch_samn_ids)
        
        try:
            # Search the BioSample database
            handle = Entrez.esearch(db="biosample", term=id_string, retmax=batch_size)
            biosample_record = Entrez.read(handle)
            handle.close()
            biosample_uids = biosample_record["IdList"]

            if not biosample_uids:
                print(f"  Processed batch {i // batch_size + 1}: No valid BioSamples found.")
                continue

            # Link from BioSample to the SRA database
            handle = Entrez.elink(dbfrom="biosample", db="sra", id=",".join(biosample_uids))
            link_records = Entrez.read(handle)
            handle.close()
            
            sra_uids = []
            for record in link_records:
                if 'LinkSetDb' in record and record['LinkSetDb']:
                    for link_db in record['LinkSetDb']:
                        if link_db.get('DbTo') == 'sra':
                            for link in link_db['Link']:
                                sra_uids.append(link['Id'])

            if not sra_uids:
                print(f"  Processed batch {i // batch_size + 1}: Found BioSamples but no linked SRA records.")
                continue

            # Use efetch to get the RunInfo table (CSV)
            handle = Entrez.efetch(db="sra", id=",".join(sra_uids), rettype='runinfo', retmode='text')
            
            # Decode bytes into a string before parsing
            csv_data = handle.read().decode('utf-8').strip().splitlines()
            runinfo_reader = csv.reader(csv_data)
            
            header = next(runinfo_reader)
            try:
                run_col_index = header.index('Run')
            except ValueError:
                print(f"  Processed batch {i // batch_size + 1}: Could not find 'Run' column in RunInfo.")
                continue
                
            batch_srr_ids = [row[run_col_index] for row in runinfo_reader]
            srr_accessions.extend(batch_srr_ids)
            print(f"  Processed batch {i // batch_size + 1}: Found {len(batch_srr_ids)} SRR accessions.")
            
            time.sleep(0.5)

        except Exception as e:
            print(f"An error occurred in a batch: {e}", file=sys.stderr)
            
    return srr_accessions

def write_commands_to_file(sra_ids, filename):
    """
    Writes prefetch and fastq-dump commands to a shell script file.
    """
    if not sra_ids:
        print(f"No SRA IDs to write for {filename}.")
        return

    with open(filename, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# This script downloads and converts SRA files to FASTQ.\n\n")
        
        batch_size = 500
        for i in range(0, len(sra_ids), batch_size):
            batch_sra_ids = sra_ids[i:i + batch_size]
            f.write(f"echo \"Prefetching batch {i // batch_size + 1}...\"\n")
            f.write(f"prefetch --max-size u {' '.join(batch_sra_ids)}\n\n")

        f.write("echo \"\nConverting .sra files to .fastq.gz...\"\n")
        f.write(f"fastq-dump --split-files --gzip --skip-technical --clip {' '.join(sra_ids)}\n")
    
    print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

def main():
    gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)]
    s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

    print("--- Fetching SRA accessions for gyrB data ---")
    gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
    print(f"Found a total of {len(gyrB_sra_ids)} gyrB SRA accessions.\n")
    
    print("--- Fetching SRA accessions for 16S data ---")
    s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
    print(f"Found a total of {len(s16S_sra_ids)} 16S SRA accessions.\n")
    
    write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
    write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

if __name__ == "__main__":
    main()
