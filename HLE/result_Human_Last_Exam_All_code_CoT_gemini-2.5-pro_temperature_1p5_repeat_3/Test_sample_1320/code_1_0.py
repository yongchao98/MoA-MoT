import sys

def identify_foldamer_helix():
    """
    Identifies the likely helix type for an alternating alpha/epsilon hybrid peptide
    based on established scientific literature.
    """
    
    # Step 1: Define the characteristics of the input peptide.
    peptide_class = "alternating alpha/epsilon-amino acid peptidomimetic"
    
    print(f"Analyzing the peptide: {peptide_class}")
    print("This class of molecule is a 'foldamer', known to form well-defined secondary structures.")
    print("-" * 50)
    
    # Step 2: Access a knowledge base of known foldamer structures.
    # In reality, this information comes from scientific journals. For this script,
    # we simulate this with a dictionary lookup.
    foldamer_helix_knowledge_base = {
        "alpha/beta": "11/9-helix",
        "beta": "12-helix or 14-helix",
        "alpha/gamma": "13/15-helix",
        # Key finding for long-chain hybrid peptides from literature:
        "alpha/delta_or_epsilon": "18/16-helix"
    }
    
    # Identify the correct entry in our knowledge base.
    lookup_key = "alpha/delta_or_epsilon"
    helix_type = foldamer_helix_knowledge_base[lookup_key]
    
    print("Consulting knowledge base for this peptide class...")
    print(f"Finding: Hybrid peptides with alternating alpha- and long-chain amino acids (like epsilon) form a specific, stable helix.")
    print("-" * 50)

    # Step 3: Output the result and its components as requested.
    print("Conclusion: The most likely structure is an 18/16-helix.")
    
    try:
        # Extract the numbers from the helix type string.
        number1, second_part = helix_type.split('/')
        number2 = second_part.split('-')[0]

        print("\nThe predicted helix type is described by the numbers:")
        print(f"First number: {number1}")
        print(f"Second number: {number2}")
        print(f"\nFinal Answer: The identified structure is the {number1}/{number2} helix.")
        
    except ValueError:
        print("\nCould not parse the helix type to extract numbers.")
        sys.exit(1)


if __name__ == "__main__":
    identify_foldamer_helix()