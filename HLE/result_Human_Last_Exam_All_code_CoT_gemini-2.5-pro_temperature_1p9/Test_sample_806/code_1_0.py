import textwrap

def explain_biochemical_reaction():
    """
    Explains the biochemical cascade leading to drug-induced skin blisters (SJS/TEN).
    """
    
    print("The patient's symptoms, particularly the development of skin blisters after starting a new anti-epileptic drug, strongly suggest a severe cutaneous adverse reaction (SCAR) like Stevens-Johnson Syndrome (SJS) or Toxic Epidermal Necrolysis (TEN).")
    print("-" * 70)
    print("The process is initiated by a very specific biochemical reaction:\n")

    # Step 1: The Initial Reaction
    step1_title = "Step 1: The Initiating Molecular Interaction"
    step1_desc = """
    The specific biochemical reaction that initiated the process is the **non-covalent binding of the causative drug (or its metabolite) to a specific Human Leukocyte Antigen (HLA) protein**. These HLA proteins are found on the surface of the body's antigen-presenting cells. This binding alters the shape of the HLA molecule, causing it to be incorrectly perceived as a foreign threat by the immune system. This concept is sometimes called the 'p-i' (pharmacological interaction with immune receptors) model.
    """
    print(f"**{step1_title}**")
    print(textwrap.fill(step1_desc.strip(), width=70))
    print("\n")

    # Step 2: Immune System Activation
    step2_title = "Step 2: T-Cell Activation"
    step2_desc = """
    This newly formed drug-HLA complex is directly recognized by receptors on specialized immune cells called cytotoxic T-lymphocytes (T-cells). This recognition activates the T-cells without the need for typical antigen processing.
    """
    print(f"**{step2_title}**")
    print(textwrap.fill(step2_desc.strip(), width=70))
    print("\n")

    # Step 3: Attack on Skin Cells
    step3_title = "Step 3: Widespread Keratinocyte Apoptosis (Cell Death)"
    step3_desc = """
    The activated T-cells multiply and travel to the skin. There, they release a cytotoxic protein called **Granulysin**. Granulysin is a key mediator that induces widespread apoptosis (programmed cell death) in keratinocytes, the primary cells of the epidermis (the outer layer of skin).
    """
    print(f"**{step3_title}**")
    print(textwrap.fill(step3_desc.strip(), width=70))
    print("\n")
    
    # Step 4: Final Outcome
    step4_title = "Step 4: Blister Formation"
    step4_desc = """
    The mass death of keratinocytes causes the epidermis to separate from the underlying dermis, leading to the formation of the characteristic painful blisters and sloughing of the skin.
    """
    print(f"**{step4_title}**")
    print(textwrap.fill(step4_desc.strip(), width=70))


if __name__ == "__main__":
    explain_biochemical_reaction()
