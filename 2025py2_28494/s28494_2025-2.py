#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever with Sequence Length Filtering, CSV Report, and Visualization
"""

from Bio import Entrez
from Bio import SeqIO
import csv
import os
import matplotlib.pyplot as plt  # For plotting


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid, min_length=None, max_length=None):
        """Search for all records matching a taxonomic ID with length filtering."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Fetch taxonomic information
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Create search term with length filters
            search_term = f"txid{taxid}[Organism]"
            if min_length:
                search_term += f" AND {min_length}:1000000[SLEN]"  # Minimum length
            if max_length:
                search_term += f" AND 0:{max_length}[SLEN]"  # Maximum length

            # Perform search
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name} with the specified criteria.")
                return None

            print(f"Found {count} records meeting the criteria.")

            # Store search results for later fetching
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count
            return count
        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []
        try:
            # Limit to prevent server overload
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            return list(SeqIO.parse(handle, "genbank"))  # Parse GenBank records into a list
        except Exception as e:
            print(f"Error fetching records: {e}")
            return []

    def generate_csv_report(self, records, output_file):
        """Generate a CSV report with record details."""
        try:
            with open(output_file, mode="w", newline="", encoding="utf-8") as csvfile:
                writer = csv.writer(csvfile)
                # Write CSV header
                writer.writerow(["Accession Number", "Sequence Length", "Description"])

                # Write record details
                for record in records:
                    accession = record.id
                    length = len(record.seq)
                    description = record.description
                    writer.writerow([accession, length, description])

            print(f"CSV report saved to {output_file}")
        except Exception as e:
            print(f"Error generating CSV report: {e}")

    def generate_plot(self, records, output_file):
        """Generate and save a line plot showing sequence lengths."""
        try:
            # Sort records by sequence length in descending order
            sorted_records = sorted(records, key=lambda rec: len(rec.seq), reverse=True)

            # Extract accession numbers and lengths
            accession_numbers = [record.id for record in sorted_records]
            sequence_lengths = [len(record.seq) for record in sorted_records]

            # Create the plot
            plt.figure(figsize=(10, 6))
            plt.plot(accession_numbers, sequence_lengths, marker='o', linestyle='-', color='b', label='Sequence Length')
            plt.xlabel("Accession Number")
            plt.ylabel("Sequence Length")
            plt.title("GenBank Records Sorted by Sequence Length")
            plt.xticks(rotation=90, fontsize=8)
            plt.tight_layout()  # Adjust layout to prevent clipping
            plt.legend()

            # Save the plot as a PNG file
            plt.savefig(output_file)
            print(f"Plot saved to {output_file}")
            plt.close()
        except Exception as e:
            print(f"Error generating plot: {e}")


def main():
    # Get user credentials for NCBI
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Create retriever object
    retriever = NCBIRetriever(email, api_key)

    # Get taxonomic ID and length filters from the user
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    try:
        min_length = int(input("Enter minimum sequence length (or leave blank for no minimum): ") or 0)
        max_length = int(input("Enter maximum sequence length (or leave blank for no maximum): ") or 1000000)
    except ValueError:
        print("Invalid length provided. Using default values.")
        min_length, max_length = None, None

    # Search for records
    count = retriever.search_taxid(taxid, min_length=min_length, max_length=max_length)
    if not count:
        print("No records found. Exiting.")
        return

    # Fetch records
    print("\nFetching records...")
    records = retriever.fetch_records(start=0, max_records=10)
    if not records:
        print("No records fetched. Exiting.")
        return

    # Generate CSV report
    output_csv_file = f"taxid_{taxid}_report.csv"
    retriever.generate_csv_report(records, output_csv_file)

    # Generate Plot
    output_plot_file = f"taxid_{taxid}_plot.png"
    retriever.generate_plot(records, output_plot_file)


if __name__ == "__main__":
    main()
