import textwrap

def identify_protein_domain():
    """
    Analyzes characteristics of protein domains to identify the one
    matching the provided conservation pattern. The image shows a pattern
    of short, tandemly repeated motifs with a high density of conserved residues.
    """

    # Step 1: Define the characteristics of the pattern seen in the image.
    image_pattern_features = {
        "repeat": "Tandem (side-by-side) repeats of a short motif.",
        "length": "Each repeat is short (estimated ~20-30 residues).",
        "conservation": "High density of conserved residues (red bars) within each repeat."
    }

    # Step 2: Define the characteristics of each candidate protein domain.
    domain_characteristics = {
        "A. SH3 domain": "A single, non-repeating domain of ~60 amino acids. Recognizes proline-rich motifs.",
        "B. Leucine zipper motif": "A repeating pattern of 7 amino acids (heptad repeat) with conserved hydrophobic residues typically at positions 1 and 4. This is a sparse conservation pattern (2 in 7).",
        "C. EGF-like domain": "A single domain of ~40 amino acids defined by a pattern of six conserved Cysteine residues.",
        "D. PH domain": "A single, non-repeating domain of ~120 amino acids that binds lipids.",
        "E. WD40 repeat": "A repeat of ~40 amino acids, longer than the pattern shown, with conservation mainly at the beginning and end of the repeat.",
        "F. PDZ domain": "A single, non-repeating domain of ~90 amino acids that binds C-terminal peptides.",
        "G. Homeobox domain": "A single, non-repeating DNA-binding domain of ~60 amino acids.",
        "H. Zinc finger domain": "A small motif (~23-30 amino acids) that is frequently found in tandem arrays. Each repeat has several conserved residues (e.g., Cys2His2 for zinc binding, plus structural residues), matching the pattern of short repeats with high conservation density."
    }

    # Step 3: Find the best match by comparing image features to domain characteristics.
    print("Analyzing the visual pattern against known protein domain features:\n")
    print(f"Image Pattern: {image_pattern_features['repeat']} {image_pattern_features['conservation']}\n")
    print("------------------------------------------------------------------")

    best_match_key = None
    # The logic is that only Zinc fingers match the "tandem short repeat" and "high density" pattern.
    for key, desc in domain_characteristics.items():
        is_match = "No"
        # Heuristic check for the best match
        if "tandem arrays" in desc or ("repeating pattern" in desc and "short" in image_pattern_features['length'].lower() and "high conservation" in desc):
             is_match = "Yes"
             best_match_key = key
        
        print(f"Checking {key}:")
        print(textwrap.fill(desc, width=70))
        # Logic to identify the match based on keywords
        if "tandem" in desc and "high conservation" in desc:
            print("--> Match Found: This domain occurs in tandem repeats with high conservation, fitting the image.\n")
            best_match_key = key
        else:
            print("--> No Match: This does not fit the tandem, high-density conservation pattern.\n")


    # Step 4: Print the final conclusion
    final_answer = best_match_key[0]
    print("==================================================================")
    print("Conclusion: The pattern of short, tandemly repeated motifs with many conserved")
    print("residues is a classic signature of a C2H2 Zinc finger protein, where multiple")
    print("finger domains are arranged consecutively to bind DNA.")
    print(f"The correct option is H.")
    print("==================================================================")
    
    return final_answer

if __name__ == '__main__':
    answer = identify_protein_domain()
    print(f"<<<{answer}>>>")