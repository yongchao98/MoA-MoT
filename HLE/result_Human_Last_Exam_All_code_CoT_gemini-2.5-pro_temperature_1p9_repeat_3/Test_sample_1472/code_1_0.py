import textwrap

def find_point_group_of_complex():
    """
    Prints a step-by-step determination of the point group for
    bis(2,5-dithiahexane)copper.
    """
    
    plan_header = "### Plan to Determine the Symmetry Point Group ###"
    step1_plan = "1. Deconstruct the molecule's name: 'bis(2,5-dithiahexane)copper'."
    step2_plan = "2. Determine the likely oxidation state and coordination geometry."
    step3_plan = "3. Analyze the structure of the most stable isomer."
    step4_plan = "4. Identify all present symmetry elements to assign the point group."
    
    analysis_header = "### Step-by-Step Analysis ###"
    
    step1_analysis_title = "Step 1: The Molecule's Composition"
    step1_analysis_text = ("- Ligand: The ligand is 2,5-dithiahexane. Its structure is CH3-S-CH2-CH2-S-CH3. "
                           "It is a neutral, bidentate ligand that coordinates through its two sulfur atoms, forming a five-membered chelate ring with the metal.\n"
                           "- Metal: The central atom is Copper (Cu).\n"
                           "- Complex: The prefix 'bis' indicates two ligands are coordinated to the copper center. The formula is [Cu(2,5-dithiahexane)2].")

    step2_analysis_title = "Step 2: Copper's Oxidation State and Coordination Geometry"
    step2_analysis_text = ("- Oxidation State: The name does not specify a charge, so we assume the most common oxidation state for copper in coordination complexes, which is Cu(II). The complex is therefore the cation [Cu(dth)2]2+.\n"
                           "- Geometry: Cu(II) has a d9 electron configuration. For a 4-coordinate complex, d9 metal ions strongly favor a square planar geometry due to the Jahn-Teller effect, which removes orbital degeneracy and lowers the energy.")

    step3_analysis_title = "Step 3: Structure of the Most Stable Isomer"
    step3_analysis_text = ("- Core Structure: The complex features a central Cu atom bonded to four S atoms from the two ligands, forming a square planar CuS4 core.\n"
                           "- Chelate Rings: Each ligand's backbone (Cu-S-C-C-S) is puckered, not flat. This puckering is chiral and can be described as either lambda (λ) or delta (δ).\n"
                           "- Isomerism: With two such rings, the most stable isomer is generally the 'meso' form where the rings have opposite puckering (one λ, one δ). This arrangement minimizes steric repulsion between the ligands.")
                           
    step4_analysis_title = "Step 4: Symmetry Element Analysis"
    step4_analysis_text = ("We analyze the symmetry of this stable 'λδ' isomer:\n"
                           "- E (Identity): Always present.\n"
                           "- i (Center of Inversion): A center of inversion exists at the Cu atom. Inverting any atom through this point maps it onto an identical atom. The λ-puckered ring is transformed into the δ-puckered ring.\n"
                           "- Cn (Proper Rotation Axis): There are no C2, C3, or C4 axes. A C2 axis perpendicular to the CuS4 plane is absent because the methyl (CH3) groups on the sulfur atoms break this symmetry.\n"
                           "- σ (Mirror Plane): There are no mirror planes. The out-of-plane CH3 groups and the puckered backbones prevent any part of the molecule from being a mirror plane.\n"
                           "- Sn (Improper Rotation Axis): An inversion center 'i' is equivalent to an S2 axis. No other Sn axes are present.")

    conclusion_header = "### Conclusion ###"
    conclusion_text = "The only symmetry elements present in the molecule are the identity element (E) and a center of inversion (i). A point group defined by these two elements is Ci."
    
    final_result_text = "\nThe point group of bis(2,5-dithiahexane)copper is therefore determined to be Ci."

    # Print the entire explanation
    print(plan_header)
    print(textwrap.fill(step1_plan, width=90))
    print(textwrap.fill(step2_plan, width=90))
    print(textwrap.fill(step3_plan, width=90))
    print(textwrap.fill(step4_plan, width=90))
    print("\n" + analysis_header)
    
    print("\n" + step1_analysis_title)
    print(textwrap.fill(step1_analysis_text, width=90))
    
    print("\n" + step2_analysis_title)
    print(textwrap.fill(step2_analysis_text, width=90))
    
    print("\n" + step3_analysis_title)
    print(textwrap.fill(step3_analysis_text, width=90))

    print("\n" + step4_analysis_title)
    print(textwrap.fill(step4_analysis_text, width=90))
    
    print("\n" + conclusion_header)
    print(textwrap.fill(conclusion_text, width=90))

    print(final_result_text)

# Execute the function to print the analysis
find_point_group_of_complex()