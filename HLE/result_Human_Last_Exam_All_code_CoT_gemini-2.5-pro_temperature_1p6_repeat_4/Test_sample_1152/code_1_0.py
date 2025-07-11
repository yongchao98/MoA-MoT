def analyze_bee_data():
    """
    Analyzes experimental data on honeybees and fungi to determine the most accurate statement.
    """
    print("Analyzing the data to verify the statement: 'Fungus A and B are pathogens. Fungus C is a commensal.'\n")

    # --- Data from Experiments ---
    baseline_mortality = 10  # Mortality rate of non-infected bees in %
    
    # Experiment 1 & 3: Mortality with Fungus A and B
    mortality_fungus_A = 35  # General mortality rate with Fungus A in %
    mortality_fungus_B = 20  # General mortality rate with Fungus B in %
    
    # Experiment 4: Mortality with Fungus C
    mortality_fungus_C = 10  # General mortality rate with Fungus C in %

    # Experiment 5: Productivity with Fungus C (Buck pollen example)
    productivity_buck_not_infected_C = 45
    productivity_buck_infected_C = 60

    # --- Step 1: Evaluate Pathogenicity of Fungus A ---
    print("Step 1: Is Fungus A a pathogen?")
    print(f"Baseline mortality rate: {baseline_mortality}%")
    print(f"Mortality rate with Fungus A: {mortality_fungus_A}%")
    if mortality_fungus_A > baseline_mortality:
        print(f"Finding: {mortality_fungus_A} > {baseline_mortality}. The mortality rate is higher, so Fungus A is a pathogen.\n")
    else:
        print("Finding: The mortality rate is not higher, so Fungus A is not a pathogen.\n")

    # --- Step 2: Evaluate Pathogenicity of Fungus B ---
    print("Step 2: Is Fungus B a pathogen?")
    print(f"Baseline mortality rate: {baseline_mortality}%")
    print(f"Mortality rate with Fungus B: {mortality_fungus_B}%")
    if mortality_fungus_B > baseline_mortality:
        print(f"Finding: {mortality_fungus_B} > {baseline_mortality}. The mortality rate is higher, so Fungus B is a pathogen.\n")
    else:
        print("Finding: The mortality rate is not higher, so Fungus B is not a pathogen.\n")
        
    # --- Step 3: Evaluate Nature of Fungus C ---
    print("Step 3: Is Fungus C a commensal?")
    print("First, checking mortality rate...")
    print(f"Baseline mortality rate: {baseline_mortality}%")
    print(f"Mortality rate with Fungus C: {mortality_fungus_C}%")
    if mortality_fungus_C == baseline_mortality:
        print(f"Finding: {mortality_fungus_C} == {baseline_mortality}. The mortality rate is not increased, suggesting it is not harmful (non-pathogenic).")
        print("Next, checking productivity effect (e.g., with Buck pollen)...")
        print(f"Egg production (not infected): {productivity_buck_not_infected_C}")
        print(f"Egg production (infected with C): {productivity_buck_infected_C}")
        if productivity_buck_infected_C > productivity_buck_not_infected_C:
            print(f"Finding: {productivity_buck_infected_C} > {productivity_buck_not_infected_C}. Productivity increased.")
        print("Conclusion: Fungus C does not cause harm and may even be beneficial. Classifying it as a commensal is accurate.\n")
    else:
        print("Finding: The mortality rate changed, so it is not a commensal.\n")

analyze_bee_data()
<<<I>>>