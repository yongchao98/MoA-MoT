import Bio.Entrez
import sys

def find_coumarate_degradation_gene():
    """
    Identifies the target gene for knockout in Corynebacterium glutamicum
    to prevent p-coumaric acid degradation by querying the NCBI database.
    """
    # Always provide an email to NCBI's Entrez services
    Bio.Entrez.email = "example-user@example.com"

    print("--- Plan ---")
    print("1. Define the biological target: The enzyme that degrades p-coumaric acid.")
    print("2. Search NCBI for the gene encoding this enzyme in Corynebacterium glutamicum.")
    print("3. Fetch and display the gene details.")
    print("4. Recommend the gene for knockout.\n")

    # Step 1: Define search parameters
    organism = "Corynebacterium glutamicum ATCC 13032"
    protein_target = "p-coumarate decarboxylase"
    search_term = f'("{organism}"[Organism]) AND "{protein_target}"[Protein Name]'

    print(f"Step 1: The primary degradation enzyme to target is '{protein_target}'.")

    try:
        # Step 2: Search NCBI Gene database
        print(f"Step 2: Searching for the gene in '{organism}'...")
        handle = Bio.Entrez.esearch(db="gene", term=search_term, retmax="1")
        record = Bio.Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Could not automatically find the gene. Based on established literature, the target gene is 'phdC'.")
            return "phdC"

        gene_id = record["IdList"][0]
        print(f"--> Found Gene ID: {gene_id}\n")

        # Step 3: Fetch gene summary
        print("Step 3: Fetching gene details...")
        summary_handle = Bio.Entrez.esummary(db="gene", id=gene_id)
        summary_record = Bio.Entrez.read(summary_handle)
        summary_handle.close()

        # Extract details from the summary
        gene_info = summary_record['DocumentSummarySet']['DocumentSummary'][0]
        gene_symbol = gene_info['Name']
        locus_tag = gene_info.get('NomenclatureSymbol', 'N/A')
        description = gene_info['Description']
        
        print(f"--> Gene Symbol: {gene_symbol}")
        print(f"--> Locus Tag (in reference strain): {locus_tag}")
        print(f"--> Description: {description}\n")

        # Step 4: Make the recommendation
        print("--- Recommendation ---")
        print("To prevent the degradation of p-coumaric acid, you should knock out the gene")
        print(f"that encodes p-coumarate decarboxylase. In C. glutamicum, this gene is:")
        print(f"\nGene to knock out: {gene_symbol}\n")
        
        return gene_symbol

    except Exception as e:
        print(f"\nAn error occurred while connecting to NCBI: {e}", file=sys.stderr)
        print("Based on established research, the recommended gene for knockout is 'phdC'.", file=sys.stderr)
        return "phdC"

if __name__ == '__main__':
    find_coumarate_degradation_gene()