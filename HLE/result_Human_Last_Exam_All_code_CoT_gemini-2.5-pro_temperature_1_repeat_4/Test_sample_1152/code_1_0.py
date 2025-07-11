def analyze_fungi_effects():
    """
    Analyzes experimental data to classify fungi based on their effect on honeybee mortality.
    """
    # --- Data from the experiments ---
    baseline_mortality = 10  # Mortality rate of non-infected honeybees in %
    fungus_a_mortality = 35  # General mortality rate for Fungus A infection in %
    fungus_b_mortality = 20  # General mortality rate for Fungus B infection in %
    fungus_c_mortality = 10  # General mortality rate for Fungus C infection in %
    
    # Productivity data for Fungus C with Buck pollen
    fungus_c_eggs_uninfected = 45
    fungus_c_eggs_infected = 60

    print("--- Fungi Analysis ---")
    print(f"Baseline mortality rate for healthy bees: {baseline_mortality}%")
    print("-" * 30)

    # --- Analysis for Fungus A ---
    print("Analyzing Fungus A:")
    print(f"Comparing mortality with Fungus A ({fungus_a_mortality}%) to baseline ({baseline_mortality}%).")
    if fungus_a_mortality > baseline_mortality:
        print(f"Result: {fungus_a_mortality} > {baseline_mortality}. The mortality rate increases, therefore Fungus A is a pathogen.")
    else:
        print(f"Result: {fungus_a_mortality} <= {baseline_mortality}. The mortality rate does not increase.")
    print("-" * 30)

    # --- Analysis for Fungus B ---
    print("Analyzing Fungus B:")
    print(f"Comparing mortality with Fungus B ({fungus_b_mortality}%) to baseline ({baseline_mortality}%).")
    if fungus_b_mortality > baseline_mortality:
        print(f"Result: {fungus_b_mortality} > {baseline_mortality}. The mortality rate increases, therefore Fungus B is a pathogen.")
    else:
        print(f"Result: {fungus_b_mortality} <= {baseline_mortality}. The mortality rate does not increase.")
    print("-" * 30)

    # --- Analysis for Fungus C ---
    print("Analyzing Fungus C:")
    print(f"Comparing mortality with Fungus C ({fungus_c_mortality}%) to baseline ({baseline_mortality}%).")
    if fungus_c_mortality > baseline_mortality:
        print(f"Result: {fungus_c_mortality} > {baseline_mortality}. The mortality rate increases, therefore Fungus C is a pathogen.")
    else:
        print(f"Result: {fungus_c_mortality} == {baseline_mortality}. The mortality rate does not increase.")
        print("Analyzing productivity for Fungus C...")
        print(f"Egg count for uninfected bees (Buck pollen): {fungus_c_eggs_uninfected}")
        print(f"Egg count for bees infected with Fungus C (Buck pollen): {fungus_c_eggs_infected}")
        if fungus_c_eggs_infected > fungus_c_eggs_uninfected:
            print(f"Result: {fungus_c_eggs_infected} > {fungus_c_eggs_uninfected}. Productivity increases.")
        print("Since it does not increase mortality and may increase productivity, Fungus C is a commensal organism.")
    print("-" * 30)
    
    print("\nFinal Conclusion: Based on the analysis, Fungus A and B are pathogens, while Fungus C is a commensal.")

# Run the analysis
analyze_fungi_effects()
<<<I>>>