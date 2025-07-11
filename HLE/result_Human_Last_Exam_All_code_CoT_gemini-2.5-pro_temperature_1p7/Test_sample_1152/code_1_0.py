def analyze_bee_experiments():
    """
    Analyzes the provided honeybee experiment data to determine the correct statement.
    """

    # Data from the experiments
    control_mortality = 10  # Mortality rate of non-infected bees

    # Experiment 1: Fungus A mortality
    # We take the general case, not the special sunflower case, to determine pathogenicity.
    fungus_A_mortality = 35

    # Experiment 3: Fungus B mortality
    fungus_B_mortality = 20

    # Experiment 4: Fungus C mortality
    fungus_C_mortality = 10
    
    # Experiment 5: Fungus C productivity (examples)
    productivity_not_infected = {'buck': 45, 'lavender': 30}
    productivity_infected_C = {'buck': 60, 'lavender': 50}

    # Step 1: Analyze Fungus A
    print("--- Analysis of Fungus A ---")
    is_pathogen_A = fungus_A_mortality > control_mortality
    print(f"Comparing Fungus A mortality rate ({fungus_A_mortality}%) with non-infected mortality rate ({control_mortality}%).")
    if is_pathogen_A:
        print("Result: Fungus A mortality rate is higher, therefore it is a pathogen.")
    else:
        print("Result: Fungus A mortality rate is not higher, therefore it is not a pathogen.")
    print("-" * 20)

    # Step 2: Analyze Fungus B
    print("--- Analysis of Fungus B ---")
    is_pathogen_B = fungus_B_mortality > control_mortality
    print(f"Comparing Fungus B mortality rate ({fungus_B_mortality}%) with non-infected mortality rate ({control_mortality}%).")
    if is_pathogen_B:
        print("Result: Fungus B mortality rate is higher, therefore it is a pathogen.")
    else:
        print("Result: Fungus B mortality rate is not higher, therefore it is not a pathogen.")
    print("-" * 20)

    # Step 3: Analyze Fungus C
    print("--- Analysis of Fungus C ---")
    is_pathogen_C = fungus_C_mortality > control_mortality
    print(f"Comparing Fungus C mortality rate ({fungus_C_mortality}%) with non-infected mortality rate ({control_mortality}%).")
    if is_pathogen_C:
        print("Result: Fungus C mortality rate is higher, therefore it is a pathogen.")
    else:
        print("Result: Fungus C does not increase mortality. Let's check productivity.")
        # Check productivity effect
        if productivity_infected_C['buck'] > productivity_not_infected['buck']:
            print(f"Productivity with buck pollen increased from {productivity_not_infected['buck']} to {productivity_infected_C['buck']} eggs.")
        if productivity_infected_C['lavender'] > productivity_not_infected['lavender']:
            print(f"Productivity with lavender pollen increased from {productivity_not_infected['lavender']} to {productivity_infected_C['lavender']} eggs.")
        print("Result: Fungus C does not cause harm and may be beneficial. It can be classified as a commensal or mutualist.")
    print("-" * 20)
    
    # Step 4: Final Conclusion based on analysis
    print("\n--- Final Conclusion ---")
    print("Our analysis shows:")
    print("1. Fungus A is a pathogen.")
    print("2. Fungus B is a pathogen.")
    print("3. Fungus C is a commensal (or mutualist).")
    print("\nThis corresponds to option I.")

analyze_bee_experiments()

print("\n<<<I>>>")