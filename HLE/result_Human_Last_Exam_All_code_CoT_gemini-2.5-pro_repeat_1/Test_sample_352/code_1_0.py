def design_mutagenesis():
    """
    Designs and explains a site-directed mutagenesis experiment to
    relieve an inhibitory, negatively charged amino acid patch.
    """
    # Define the original amino acids and their positions in the protein.
    original_patch = {
        47: {"code": "S", "name": "Serine"},
        48: {"code": "E", "name": "Glutamate"},
        49: {"code": "E", "name": "Glutamate"},
        50: {"code": "D", "name": "Aspartate"}
    }

    # Define the proposed replacement amino acid.
    # Alanine (A) is chosen to remove negative charge and prevent phosphorylation.
    replacement_aa = {"code": "A", "name": "Alanine"}

    print("Experiment Plan: Site-Directed Mutagenesis of Protein X")
    print("-" * 55)
    print("Objective: To eliminate the autoinhibitory effect of the negatively")
    print("           charged S-E-E-D patch at positions 47-50.\n")
    print("Strategy: Replace all four residues with Alanine (A) to create a")
    print("          neutral, non-phosphorylatable, and sterically minimal patch.\n")

    print("Proposed mutations:")
    
    original_sequence = ""
    replacement_sequence = ""
    
    # Iterate through the sorted positions and print the mutation details.
    for position in sorted(original_patch.keys()):
        original_code = original_patch[position]["code"]
        original_name = original_patch[position]["name"]
        
        # Build the sequences for the final output line
        original_sequence += original_code
        replacement_sequence += replacement_aa["code"]

        print(f"  Position {position}: Replace {original_code} ({original_name}) with {replacement_aa['code']} ({replacement_aa['name']})")

    print("\nSummary of change:")
    print(f"The native sequence segment '{original_sequence}' at positions 47-50 will be mutated to '{replacement_sequence}'.")

# Execute the function to print the plan.
design_mutagenesis()