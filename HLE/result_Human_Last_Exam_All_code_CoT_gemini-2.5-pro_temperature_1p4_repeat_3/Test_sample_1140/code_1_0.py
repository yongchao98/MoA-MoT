import textwrap

def explain_peptide_synthesis_choice():
    """
    Analyzes and explains the best method for synthesizing a specific peptide.
    """
    # 1. Define the problem
    peptide_length = 100  # amino acids
    unnatural_aa = "Azido phenylalanine (AzF)"
    sequence_snippet = "...KAVCLXVIGATR..."

    print("--- Analysis of Peptide Synthesis Problem ---")
    print(f"Target: A {peptide_length} amino acid peptide.")
    print(f"Special Feature: Contains the unnatural amino acid (UAA) '{unnatural_aa}' represented by 'X'.")
    print(f"Sequence Information: The middle of the sequence is '{sequence_snippet}', which includes a Cysteine (C).\n")

    # 2. Evaluate options
    print("--- Evaluation of Synthesis Techniques ---")
    
    # Solid-Phase Peptide Synthesis (SPPS)
    print("\n1. Solid-Phase Peptide Synthesis (SPPS):")
    spps_desc = """
    This is the standard method for chemical peptide synthesis. While it can incorporate UAAs like AzF, its efficiency drops significantly for peptides longer than ~50-60 amino acids. Synthesizing a 100aa peptide directly would result in very low yields and a difficult purification process.
    Conclusion: Not suitable for the full-length synthesis.
    """
    print(textwrap.indent(textwrap.dedent(spps_desc).strip(), "   "))

    # Genetic Code Expansion
    print("\n2. Genetic Code Expansion:")
    gce_desc = """
    This biological technique engineers an organism (like E. coli) to incorporate a UAA during protein translation. It is powerful for producing large quantities of UAA-containing proteins. However, it requires significant molecular biology expertise to develop and optimize the orthogonal tRNA/synthetase system.
    Conclusion: A viable but complex biological approach.
    """
    print(textwrap.indent(textwrap.dedent(gce_desc).strip(), "   "))

    # Native Chemical Ligation (NCL)
    print("\n3. Native Chemical Ligation (NCL):")
    ncl_desc = f"""
    NCL is a chemical method used to join two smaller, unprotected peptide fragments. One fragment needs a C-terminal thioester, and the other needs an N-terminal Cysteine. This is the key insight for this problem. The 100aa peptide can be split into two ~50aa fragments at the Cysteine residue found in the sequence '{sequence_snippet}'.
    - Fragment 1 (...KAV) can be synthesized as a C-terminal thioester.
    - Fragment 2 (CLXVIGATR...) has a natural N-terminal Cysteine and contains the UAA 'X'.
    Both ~50aa fragments can be efficiently made using SPPS, incorporating the UAA in the process. The fragments are then ligated to form the full, pure 100aa peptide.
    Conclusion: This is the most practical and efficient chemical synthesis strategy.
    """
    print(textwrap.indent(textwrap.dedent(ncl_desc).strip(), "   "))

    # 3. Final Conclusion and Equation
    print("\n--- Final Recommendation ---")
    print("Native Chemical Ligation (NCL) is the most helpful technique for this synthesis.")
    print("It leverages the strengths of SPPS for creating manageable fragments (including the UAA) and overcomes the length limitations of direct SPPS.\n")
    
    print("The synthesis can be represented by the following chemical equation:")
    fragment_1 = "(Peptide_1)...KAV-(Thioester)"
    fragment_2 = "H2N-Cys-LXVIGATR...(Peptide_2)"
    product = "(Peptide_1)...KAVCLXVIGATR...(Peptide_2)"
    # Printing each "number" (in this case, component) of the final equation
    print("Component 1: ", fragment_1)
    print("Component 2: ", fragment_2)
    print("Product: ", product)
    print(f"\nFinal Equation: {fragment_1} + {fragment_2} ---> {product}")


if __name__ == "__main__":
    explain_peptide_synthesis_choice()