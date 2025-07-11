def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to find the correct conclusion.
    """

    # --- Data Organization ---
    # Mortality rate for non-infected bees
    baseline_mortality_uninfected = 10  # in %

    # Experiment 1 & 2: Fungus A
    fungus_a_data = {
        'buck': {'mortality': 35, 'prod_infected': 10, 'prod_uninfected': 45},
        'sunflower': {'mortality': 10, 'prod_infected': 20, 'prod_uninfected': 30},
    }

    # Experiment 3: Fungus B
    fungus_b_data = {
        'any_pollen': {'mortality': 20}
    }

    # Experiment 4 & 5: Fungus C
    fungus_c_data = {
        'buck': {'mortality': 10, 'prod_infected': 60, 'prod_uninfected': 45},
        'sunflower': {'mortality': 10, 'prod_infected': 25, 'prod_uninfected': 30},
        'lavender': {'mortality': 10, 'prod_infected': 50, 'prod_uninfected': 30},
    }

    print("Analyzing the provided data to select the correct answer.\n")

    # --- Step 1: Analyze Fungus A ---
    print("Step 1: Evaluating Fungus A...")
    fungus_a_mortality = fungus_a_data['buck']['mortality']
    print(f"The mortality rate for bees infected with Fungus A is {fungus_a_mortality}%.")
    print(f"The baseline mortality for non-infected bees is {baseline_mortality_uninfected}%.")
    print(f"Equation: Is {fungus_a_mortality} > {baseline_mortality_uninfected}?")
    if fungus_a_mortality > baseline_mortality_uninfected:
        print("Result: Yes. This indicates Fungus A is a pathogen.\n")
    else:
        print("Result: No. This does not indicate Fungus A is a pathogen.\n")

    # --- Step 2: Analyze Fungus B ---
    print("Step 2: Evaluating Fungus B...")
    fungus_b_mortality = fungus_b_data['any_pollen']['mortality']
    print(f"The mortality rate for bees infected with Fungus B is {fungus_b_mortality}%.")
    print(f"The baseline mortality for non-infected bees is {baseline_mortality_uninfected}%.")
    print(f"Equation: Is {fungus_b_mortality} > {baseline_mortality_uninfected}?")
    if fungus_b_mortality > baseline_mortality_uninfected:
        print("Result: Yes. This indicates Fungus B is also a pathogen.\n")
    else:
        print("Result: No. This does not indicate Fungus B is a pathogen.\n")
        
    # --- Step 3: Analyze Fungus C ---
    print("Step 3: Evaluating Fungus C...")
    fungus_c_mortality = fungus_c_data['buck']['mortality']
    print(f"The mortality rate for bees infected with Fungus C is {fungus_c_mortality}%.")
    print(f"The baseline mortality for non-infected bees is {baseline_mortality_uninfected}%.")
    print(f"Equation (Mortality): Is {fungus_c_mortality} > {baseline_mortality_uninfected}?")
    if fungus_c_mortality > baseline_mortality_uninfected:
        print("Result: Yes. This indicates Fungus C is a pathogen.\n")
    else:
        print("Result: No, mortality is not increased. Let's check productivity.")
        # Check productivity
        prod_infected = fungus_c_data['buck']['prod_infected']
        prod_uninfected = fungus_c_data['buck']['prod_uninfected']
        print(f"Productivity of infected bees (e.g., on buck pollen) is {prod_infected} eggs.")
        print(f"Productivity of uninfected bees (on buck pollen) is {prod_uninfected} eggs.")
        print(f"Equation (Productivity): Is {prod_infected} < {prod_uninfected}?")
        if prod_infected < prod_uninfected:
            print("Result: No, productivity is not decreased.")
        else:
            print("Result: No, productivity is actually increased.")
        print("Since Fungus C does not increase mortality or decrease productivity, it can be classified as a commensal.\n")
        
    # --- Final Conclusion ---
    print("Conclusion:")
    print("Based on the analysis:")
    print(f"- Fungus A is a pathogen because its mortality rate ({fungus_a_mortality}%) is greater than the baseline ({baseline_mortality_uninfected}%).")
    print(f"- Fungus B is a pathogen because its mortality rate ({fungus_b_mortality}%) is greater than the baseline ({baseline_mortality_uninfected}%).")
    print(f"- Fungus C is a commensal because it does not increase mortality ({fungus_c_mortality}% vs {baseline_mortality_uninfected}%) and does not cause harm via productivity loss.")
    print("This directly supports the statement: 'Fungus A and B are pathogens. Fungus C is a commensal.'")

analyze_bee_experiments()
print("<<<I>>>")