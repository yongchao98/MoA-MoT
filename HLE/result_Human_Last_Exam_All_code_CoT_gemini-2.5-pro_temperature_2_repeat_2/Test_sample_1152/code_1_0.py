def analyze_bee_data():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    # Baseline data from Experiment 1
    baseline_mortality_rate = 10  # % for non-infected bees

    # Experiment 1 & 2 Data: Fungus A
    fungus_a_mortality_rate_general = 35 # %
    
    # Experiment 3 Data: Fungus B
    fungus_b_mortality_rate = 20 # %

    # Experiment 4 & 5 Data: Fungus C
    fungus_c_mortality_rate = 10 # %
    eggs_uninfected_buck = 45
    eggs_infected_c_buck = 60
    eggs_uninfected_other = 30
    eggs_infected_c_other = 50

    print("Step 1: Analyzing Pathogenicity of Fungi")
    print("----------------------------------------")
    
    # Analysis for Fungus A
    print("\nAnalyzing Fungus A:")
    print(f"The baseline mortality rate for non-infected bees is {baseline_mortality_rate}%.")
    print(f"With Fungus A infection, the mortality rate increases to {fungus_a_mortality_rate_general}%.")
    if fungus_a_mortality_rate_general > baseline_mortality_rate:
        print(f"Since {fungus_a_mortality_rate_general} > {baseline_mortality_rate}, Fungus A causes increased mortality and is a pathogen.")
    else:
        print("Fungus A is NOT a pathogen.")
        
    # Analysis for Fungus B
    print("\nAnalyzing Fungus B:")
    print(f"The baseline mortality rate for non-infected bees is {baseline_mortality_rate}%.")
    print(f"With Fungus B infection, the mortality rate increases to {fungus_b_mortality_rate}%.")
    if fungus_b_mortality_rate > baseline_mortality_rate:
        print(f"Since {fungus_b_mortality_rate} > {baseline_mortality_rate}, Fungus B causes increased mortality and is a pathogen.")
    else:
        print("Fungus B is NOT a pathogen.")

    # Analysis for Fungus C
    print("\nAnalyzing Fungus C:")
    print(f"The baseline mortality rate for non-infected bees is {baseline_mortality_rate}%.")
    print(f"With Fungus C infection, the mortality rate is {fungus_c_mortality_rate}%.")
    if fungus_c_mortality_rate > baseline_mortality_rate:
        print(f"Since {fungus_c_mortality_rate} is not greater than {baseline_mortality_rate}, Fungus C does not increase mortality.")
    else:
        print(f"Since {fungus_c_mortality_rate} is equal to the baseline {baseline_mortality_rate}, Fungus C does not increase mortality.")

    print("\nChecking productivity for Fungus C:")
    print(f"For bees fed on most pollens, egg count went from {eggs_uninfected_other} (not infected) to {eggs_infected_c_other} (infected with Fungus C).")
    print("This shows an increase in productivity.")
    print("Since Fungus C does not cause increased mortality or a consistent decrease in productivity, it is best described as a commensal (or even mutualist), not a pathogen.")

    print("\n--- Conclusion ---")
    print("Based on the analysis:")
    print("1. Fungus A is a pathogen.")
    print("2. Fungus B is a pathogen.")
    print("3. Fungus C is a commensal.")
    print("This directly corresponds to Answer Choice I.")

analyze_bee_data()