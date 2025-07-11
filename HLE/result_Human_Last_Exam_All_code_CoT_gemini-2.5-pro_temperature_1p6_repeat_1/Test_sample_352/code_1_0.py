import textwrap

def suggest_mutagenesis():
    """
    Analyzes the protein x inhibitory patch and suggests the best replacement amino acids.
    """
    
    # Define the original amino acid patch and its properties
    original_patch = {
        47: {'code': 'S', 'name': 'Serine', 'problem': 'A known phosphorylation site, which adds a strong negative charge.'},
        48: {'code': 'E', 'name': 'Glutamate', 'problem': 'Acidic and negatively charged at physiological pH.'},
        49: {'code': 'E', 'name': 'Glutamate', 'problem': 'Acidic and negatively charged at physiological pH.'},
        50: {'code': 'D', 'name': 'Aspartate', 'problem': 'Acidic and negatively charged at physiological pH.'}
    }
    
    # Define the proposed replacement amino acid
    replacement_aa = {'code': 'A', 'name': 'Alanine'}
    
    # --- Print the analysis and recommendation ---
    
    print("--- Analysis of the Autoinhibitory Patch (Positions 47-50) ---")
    
    original_sequence_str = " ".join([f"{pos}{details['code']}" for pos, details in original_patch.items()])
    print(f"Original Patch Sequence: {original_sequence_str}")
    print("\nThe problem is the highly negative charge of this patch, which causes autoinhibition:")
    for pos, details in original_patch.items():
        print(f"- Position {pos} ({details['name']}, {details['code']}): {details['problem']}")
        
    print("\n--- Proposed Solution: Alanine Scanning Mutagenesis ---")
    
    explanation = f"""
The best replacement for all four residues is Alanine (A). This is a standard and effective strategy for several reasons:
1. Neutral Charge: Alanine is a non-polar, uncharged amino acid. It will completely neutralize the negative charge of the patch.
2. Eliminates Phosphorylation: Replacing Serine (S47) with Alanine (A) removes the hydroxyl group necessary for phosphorylation, thus eliminating this key source of negative regulation.
3. Minimal Structural Impact: Alanine is the second smallest amino acid. Its small, inert side chain (-CH3) removes the original functionality without introducing steric bulk or new chemical interactions that might complicate the results. This is ideal for an IDR which lacks a fixed structure.
"""
    print(textwrap.dedent(explanation).strip())
    
    print("--- Final Recommendation for Mutagenesis ---")
    
    original_str = ""
    mutated_str = ""
    for pos, details in original_patch.items():
        original_str += f"{pos}{details['code']} "
        mutated_str += f"{pos}{replacement_aa['code']} "

    print(f"To relieve inhibition, mutate the original patch:")
    print(f"Original: {original_str.strip()}")
    print(f"To:")
    print(f"Mutated:  {mutated_str.strip()}")
    print("\nThis creates a neutral 'AAAA' patch, directly testing the hypothesis that the negative charge of the 'SEED' motif is inhibitory.")

# Execute the function to print the suggestion
if __name__ == "__main__":
    suggest_mutagenesis()