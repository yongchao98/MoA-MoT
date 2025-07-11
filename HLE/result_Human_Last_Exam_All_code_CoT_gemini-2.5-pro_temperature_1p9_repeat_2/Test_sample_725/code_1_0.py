def analyze_covalency():
    """
    Analyzes and explains the relative covalency between CeF₆²⁻ and CeCl₆²⁻
    based on the principle of orbital overlap.
    """
    
    # Define the core principle
    principle = "In chemical bonding, the degree of covalency is directly proportional to the extent of orbital overlap between the bonding atoms. Greater overlap signifies stronger covalency."

    # Define the complexes and orbitals from the problem
    complex_fluoride = "CeF₆²⁻"
    complex_chloride = "CeCl₆²⁻"
    metal_ion_oxidation_state = "IV"
    metal_orbital = "4f"
    fluorine_orbital = "2p"
    chlorine_orbital = "3p"

    print("Step 1: State the Chemical Principle")
    print(principle)
    print("-" * 60)

    print("Step 2: Apply the Principle to the Given Information")
    print(f"The problem states that the orbital overlap in {complex_fluoride} is greater than in {complex_chloride}.")
    print("This refers to the overlap between the central Cerium(IV) ion's "
          f"{metal_orbital} orbitals and the ligand orbitals.")
    print("-" * 60)

    print("Step 3: Draw the Final Conclusion")
    print(f"Given the greater orbital overlap, {complex_fluoride} will display STRONGER covalency compared to {complex_chloride}.")
    print("\nThis conclusion is based on the enhanced mixing of the Ce({metal_ion_oxidation_state}) {metal_orbital} orbital with the F {fluorine_orbital} orbital in the {complex_fluoride} complex, as compared to the mixing with the Cl {chlorine_orbital} orbital in the {complex_chloride} complex.")


# Run the analysis
analyze_covalency()