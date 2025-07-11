import textwrap

def find_integrin_binding_peptide():
    """
    Evaluates a list of RGD-containing peptides to identify the one most likely to
    bind an integrin receptor in an in vitro assay.
    """

    analysis = {
        'A. RGDMAA': "Not a well-documented high-affinity sequence from a major natural ligand.",
        'B. RGDSPSS': "Contains the 'RGDSP' motif. This sequence is part of the cell-binding domain of fibronectin, a major integrin ligand. Peptides with this motif (e.g., GRGDSP) are classic examples used in countless in vitro assays to demonstrate potent integrin binding.",
        'C. RGDLTTP': "Not a well-documented high-affinity sequence from a major natural ligand.",
        'D. RGDQVSK': "Not a well-documented high-affinity sequence from a major natural ligand.",
        'E. RGDARGG': "Not a well-documented high-affinity sequence from a major natural ligand."
    }

    best_choice_key = "B. RGDSPSS"
    conclusion_text = f"""
    The core binding motif for many integrin receptors is the RGD (Arginine-Glycine-Aspartate) sequence, present in all options.

    However, the flanking amino acids critically determine binding affinity and specificity. We need to identify the peptide sequence that is known from experimental data to be a potent binder.

    Let's analyze the options:
    {'='*25}
    """
    print(textwrap.dedent(conclusion_text).strip())

    for choice, reason in analysis.items():
        print(f"{choice}: {reason}")

    final_conclusion = f"""
    {'='*25}
    The peptide '{best_choice_key}' is derived from the cell-binding domain of fibronectin. Due to its origin from a natural, high-affinity ligand, it is the most well-studied and most likely peptide among the choices to bind strongly to an integrin receptor in an in vitro assay.
    """
    print(textwrap.dedent(final_conclusion).strip())

    # As requested, outputting the 'components' of the final answer.
    print("\n--- Final Answer Breakdown ---")
    print(f"Choice: B")
    print(f"Peptide Sequence: RGDSPSS")
    print(f"Reasoning: Contains the 'RGDSP' motif from the fibronectin binding site, a known high-affinity target for integrins.")

find_integrin_binding_peptide()