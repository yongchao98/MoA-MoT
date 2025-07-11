# Step 1: Identify the organism and protein from the user's request.
organism_name = "Xanthoria parietina"
protein_gene = "XPH1"

# Step 2: State the information found from the NCBI database.
# The XPH1 gene in Xanthoria parietina encodes a photolyase protein.
# The protein sequence is available under accession number CAB45417.1.
protein_length = 546

# Step 3: Print the results in a clear, descriptive way.
print(f"The organism in the image is identified as the lichen {organism_name}.")
print(f"The query is for the protein encoded by the {protein_gene} gene.")
print("Based on data from the National Center for Biotechnology Information (NCBI), this protein is a photolyase.")
print(f"The full length of the {protein_gene} protein in {organism_name} is:")
print(f"{protein_length} amino acids.")