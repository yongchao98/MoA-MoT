def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, pollen, and fungi to determine the correct conclusion.
    """
    
    # --- Data from the problem description ---
    
    # Baseline
    non_infected_mortality = 10  # in %

    # Experiment 1 & 2: Fungus A
    fungus_a_data = {
        'buck': {'mortality': 35, 'infected_eggs': 10},
        'sunflower': {'mortality': 10, 'infected_eggs': 20},
        'lavender': {'mortality': 35, 'infected_eggs': 10},
        'canola': {'mortality': 35, 'infected_eggs': 8},
        'milkweed': {'mortality': 35, 'infected_eggs': 9},
        'aster': {'mortality': 35, 'infected_eggs': 11},
        'mixed': {'mortality': 35, 'infected_eggs': 12}
    }

    # Experiment 3: Fungus B
    fungus_b_data = {
        'buck': {'mortality': 20},
        'sunflower': {'mortality': 20},
        'lavender': {'mortality': 20},
        'canola': {'mortality': 20},
        'milkweed': {'mortality': 20},
        'aster': {'mortality': 20},
        'mixed': {'mortality': 20}
    }

    # Experiment 4 & 5: Fungus C
    fungus_c_data = {
        'buck': {'mortality': 10, 'infected_eggs': 60},
        'sunflower': {'mortality': 10, 'infected_eggs': 25},
        'lavender': {'mortality': 10, 'infected_eggs': 50},
        'canola': {'mortality': 10, 'infected_eggs': 50},
        'milkweed': {'mortality': 10, 'infected_eggs': 50},
        'aster': {'mortality': 10, 'infected_eggs': 50},
        'mixed': {'mortality': 10, 'infected_eggs': 52}
    }
    
    print("Analyzing the experimental data step-by-step:\n")

    # Step 1: Analyze Fungus A and B to determine if they are pathogens
    print("Part 1: Are Fungus A and B pathogens?")
    print("A pathogen causes harm, which we can measure by increased mortality.")
    print(f"The baseline mortality for non-infected bees is {non_infected_mortality}%.")
    
    fungus_a_max_mortality = max(d['mortality'] for d in fungus_a_data.values())
    print(f"The highest mortality rate for bees infected with Fungus A is {fungus_a_max_mortality}%.")
    print(f"Since {fungus_a_max_mortality} > {non_infected_mortality}, Fungus A is a pathogen.")

    fungus_b_max_mortality = max(d['mortality'] for d in fungus_b_data.values())
    print(f"The mortality rate for bees infected with Fungus B is {fungus_b_max_mortality}%.")
    print(f"Since {fungus_b_max_mortality} > {non_infected_mortality}, Fungus B is also a pathogen.\n")

    # Step 2: Analyze Fungus C to determine its nature
    print("Part 2: What is the nature of Fungus C?")
    print("A commensal organism lives with a host without causing harm.")
    
    fungus_c_mortality = fungus_c_data['buck']['mortality'] # It's 10% for all
    print(f"The mortality rate for bees infected with Fungus C is {fungus_c_mortality}%.")
    print(f"The mortality with Fungus C ({fungus_c_mortality}%) is the same as the baseline ({non_infected_mortality}%).")
    print("This means Fungus C does not cause harm to the bees and can be classified as a commensal.\n")
    
    # Step 3: Final Conclusion
    print("Final Conclusion:")
    print("Based on the analysis, the statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is fully supported by the data.")

# Execute the analysis
analyze_bee_experiments()