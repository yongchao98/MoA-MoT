import textwrap

def find_optimal_catalyst():
    """
    Selects and describes an optimal catalyst system for dual-functionality
    olefin polymerization and polyolefin hydrogenolysis based on established
    chemical principles.
    """

    # Step 1: Define the potential components for the catalyst system.
    # The dictionaries store the component and a brief justification for its consideration.
    group_iv_metals = {
        "Ti": "Titanium: High polymerization activity, but potentially lower thermal stability.",
        "Zr": "Zirconium: A well-balanced choice offering good activity and superior stability compared to Ti.",
        "Hf": "Hafnium: Excellent stability, but often lower activity than Ti or Zr."
    }

    ligand_types = {
        "Metallocene (e.g., Cp2)": "Classic for polymerization, but may lack the robustness for C-C cleavage conditions.",
        "PNP Pincer Ligand": "Tridentate, offers exceptional thermal stability and a tunable electronic environment, leaving open sites for catalysis."
    }

    support_materials = {
        "Silica (SiO2)": "A common, high-surface-area support, but relatively inert.",
        "Sulfated Alumina (SO4/Al2O3)": "A solid acid support that can help generate and stabilize the active cationic metal center required for both reactions."
    }

    # Step 2: Select the optimal combination based on the dual-functionality requirement.
    # The logic prioritizes stability and features that promote both reactions.
    optimal_metal = "Zr"
    optimal_ligand = "PNP Pincer Ligand"
    optimal_support = "Sulfated Alumina (SO4/Al2O3)"

    # Step 3: Print the final recommended catalyst system and the rationale.
    print("=" * 80)
    print("Proposed Optimal Catalyst System for Dual-Function Polymer Upcycling")
    print("=" * 80)

    print(f"\n[Component 1] Metal Center:\n    {optimal_metal} ({group_iv_metals[optimal_metal]})")
    print(f"\n[Component 2] Ligand Framework:\n    {optimal_ligand} ({ligand_types[optimal_ligand]})")
    print(f"\n[Component 3] Support Material:\n    {optimal_support} ({support_materials[optimal_support]})")

    print("\n" + "=" * 80)
    print("Final Proposed Catalyst Formulation:")
    print(f"    Zr({optimal_ligand}) / {optimal_support}")
    print("=" * 80)

    rationale_text = """
    This proposed catalyst combines the balanced reactivity of Zirconium with a
    thermally robust PNP pincer ligand. The ligand's structure provides the
    necessary stability to withstand hydrogenolysis conditions while its tunable
    electronic properties can be optimized to balance the electrophilicity needed
    for polymerization against the requirements for C-C and H-H bond activation.
    Supporting this complex on sulfated alumina, a solid acid, is key. The support
    acts as a co-catalyst, generating a highly electrophilic 'cation-like' active
    site that is known to be effective for both olefin polymerization and alkane
    C-C bond cleavage. This single-site system is therefore a prime candidate to
    bridge the gap between polymer synthesis and degradation.
    """
    print("\nRationale for this Combination:")
    print(textwrap.fill(rationale_text.strip(), width=80))

# Execute the function to get the answer.
find_optimal_catalyst()