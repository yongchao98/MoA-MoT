def analyze_bee_data():
    """
    Analyzes the experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    # --- Data Setup ---
    baseline_mortality = 10  # %

    # Fungus A Data
    fungus_A = {
        'mortality': {'sunflower': 10, 'other': 35},
        'pathogenicity': None
    }

    # Fungus B Data
    fungus_B = {
        'mortality': {'all_pollens': 20},
        'pathogenicity': None
    }

    # Fungus C Data
    fungus_C = {
        'mortality': {'all_pollens': 10},
        'pathogenicity': None
    }

    # --- Analysis ---

    # Step 1: Determine pathogenicity of Fungus A
    # A pathogen increases mortality above the baseline.
    mortality_A = fungus_A['mortality']['other']
    if mortality_A > baseline_mortality:
        fungus_A['pathogenicity'] = 'Pathogen'
    else:
        fungus_A['pathogenicity'] = 'Not a Pathogen'
    
    print(f"Analysis for Fungus A:")
    print(f"Baseline mortality: {baseline_mortality}%")
    print(f"Mortality with Fungus A (most pollens): {mortality_A}%")
    print(f"Since {mortality_A} > {baseline_mortality}, Fungus A is a {fungus_A['pathogenicity']}.\n")

    # Step 2: Determine pathogenicity of Fungus B
    mortality_B = fungus_B['mortality']['all_pollens']
    if mortality_B > baseline_mortality:
        fungus_B['pathogenicity'] = 'Pathogen'
    else:
        fungus_B['pathogenicity'] = 'Not a Pathogen'

    print(f"Analysis for Fungus B:")
    print(f"Baseline mortality: {baseline_mortality}%")
    print(f"Mortality with Fungus B: {mortality_B}%")
    print(f"Since {mortality_B} > {baseline_mortality}, Fungus B is a {fungus_B['pathogenicity']}.\n")

    # Step 3: Determine nature of Fungus C
    # A commensal does not cause harm (i.e., does not increase mortality).
    mortality_C = fungus_C['mortality']['all_pollens']
    if mortality_C > baseline_mortality:
        fungus_C['pathogenicity'] = 'Pathogen'
    else:
        # If mortality is not increased, it's not a pathogen. "Commensal" is the best fit.
        fungus_C['pathogenicity'] = 'Commensal'

    print(f"Analysis for Fungus C:")
    print(f"Baseline mortality: {baseline_mortality}%")
    print(f"Mortality with Fungus C: {mortality_C}%")
    print(f"Since {mortality_C} is not greater than {baseline_mortality}, Fungus C is not a pathogen. It is best described as a {fungus_C['pathogenicity']}.\n")

    # Step 4: Formulate the final conclusion based on the analysis
    print("--- Final Conclusion ---")
    print(f"The analysis shows that Fungus A and Fungus B are both pathogens because they increase bee mortality.")
    print(f"Fungus C is a commensal because it does not increase mortality.")
    print("This corresponds to answer choice I.")


analyze_bee_data()