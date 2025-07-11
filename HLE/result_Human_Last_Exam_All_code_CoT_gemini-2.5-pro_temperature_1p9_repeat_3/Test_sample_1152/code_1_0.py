def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    # --- Data Storage ---
    # Baseline mortality for non-infected bees
    baseline_mortality = 10

    # Experiment 1 & 2: Fungus A (intestinal pathogen)
    fungus_a_mortality = {
        'buck': 35, 'sunflower': 10, 'lavender': 35, 'canola': 35, 
        'milkweed': 35, 'aster': 35, 'mixed': 35
    }
    
    # Experiment 3: Fungus B (surface pathogen)
    fungus_b_mortality = {
        'buck': 20, 'sunflower': 20, 'lavender': 20, 'canola': 20, 
        'milkweed': 20, 'aster': 20, 'mixed': 20
    }
    
    # Experiment 4 & 5: Fungus C (intestinal)
    fungus_c_mortality = {
        'buck': 10, 'sunflower': 10, 'lavender': 10, 'canola': 10, 
        'milkweed': 10, 'aster': 10, 'mixed': 10
    }
    fungus_c_productivity = {
        'buck': {'not_infected': 45, 'infected': 60},
        'lavender': {'not_infected': 30, 'infected': 50},
    }

    print("Step 1: Analyzing Pathogenicity of Fungus A and B")
    print("="*50)
    
    # Analyze Fungus A
    mortality_A_buck = fungus_a_mortality['buck']
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")
    print(f"With Fungus A infection (and buck pollen), the mortality rate is {mortality_A_buck}%.")
    print(f"Since {mortality_A_buck} > {baseline_mortality}, Fungus A is a pathogen.")
    
    print("-" * 20)
    
    # Analyze Fungus B
    mortality_B_buck = fungus_b_mortality['buck']
    print(f"With Fungus B infection (and buck pollen), the mortality rate is {mortality_B_buck}%.")
    print(f"Since {mortality_B_buck} > {baseline_mortality}, Fungus B is also a pathogen.")

    print("\nStep 2: Analyzing the Nature of Fungus C")
    print("="*50)

    # Analyze Fungus C mortality
    mortality_C_buck = fungus_c_mortality['buck']
    print(f"With Fungus C infection, the mortality rate is {mortality_C_buck}%.")
    print(f"This is equal to the baseline mortality ({baseline_mortality}%), indicating it does not kill the bees.")

    # Analyze Fungus C productivity
    prod_c_uninfected = fungus_c_productivity['buck']['not_infected']
    prod_c_infected = fungus_c_productivity['buck']['infected']
    print(f"\nFurther, let's look at productivity (egg count) with Fungus C.")
    print(f"A non-infected colony on buck pollen laid {prod_c_uninfected} eggs.")
    print(f"A Fungus C-infected colony on buck pollen laid {prod_c_infected} eggs.")
    print(f"Since productivity increased from {prod_c_uninfected} to {prod_c_infected}, the fungus provides a benefit instead of harm.")
    print("Conclusion: Fungus C is not a pathogen and can be described as a commensal or symbiotic organism.")

    print("\nFinal Conclusion:")
    print("Based on the analysis, 'Fungus A and B are pathogens. Fungus C is a commensal.' is the most accurate statement.")

analyze_bee_experiments()
<<<I>>>