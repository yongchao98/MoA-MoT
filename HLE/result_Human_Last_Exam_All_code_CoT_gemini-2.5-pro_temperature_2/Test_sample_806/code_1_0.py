import textwrap

def explain_reaction_pathway():
    """
    This function explains the biochemical cascade leading to medication-induced skin blisters (SJS/TEN).
    """

    # --- Define the key players based on the clinical scenario ---
    trigger_drug = "Anticonvulsant Drug (e.g., Lamotrigine)"
    patient_receptor = "Human Leukocyte Antigen (HLA) Molecule"
    immune_cell = "Cytotoxic T-lymphocyte (T-cell)"
    target_cell = "Keratinocyte (Skin Cell)"

    print("--- Step-by-Step Biochemical Explanation ---")

    # 1. The initiating reaction
    step1_text = f"The process begins when the {trigger_drug} molecule binds directly and non-covalently to a specific " \
                 f"{patient_receptor} on the surface of the patient's cells. This binding alters the shape of the HLA molecule."
    print("\n[STEP 1: The Initiating Reaction]")
    print(textwrap.fill(step1_text, 80))

    # 2. Immune system activation
    step2_text = f"This newly formed [Drug + {patient_receptor}] complex is now recognized as 'foreign' by the patient's " \
                 f"immune system, specifically by a {immune_cell}."
    print("\n[STEP 2: Immune Cell Recognition]")
    print(textwrap.fill(step2_text, 80))

    # 3. Cytotoxic attack
    step3_text = f"The recognition event activates the {immune_cell}, causing it to release powerful cytotoxic proteins " \
                 f"(like granulysin). These proteins target and induce apoptosis (programmed cell death) in the {target_cell}s."
    print("\n[STEP 3: Cytotoxic Attack]")
    print(textwrap.fill(step3_text, 80))

    # 4. Clinical manifestation
    step4_text = f"The widespread death of {target_cell}s leads to the separation of the epidermis from the dermis, " \
                 f"which manifests clinically as the severe skin blisters and sloughing seen in the patient."
    print("\n[STEP 4: Clinical Result]")
    print(textwrap.fill(step4_text, 80))

    # 5. Final Answer
    final_answer = "The specific biochemical reaction that initiated the process was the non-covalent binding " \
                   "of the drug to a specific Human Leukocyte Antigen (HLA) molecule."
    print("\n-------------------------------------------")
    print("CONCLUSION: The initiating biochemical reaction is:")
    print(textwrap.fill(final_answer, 80))
    print("-------------------------------------------")

# Execute the explanation
explain_reaction_pathway()