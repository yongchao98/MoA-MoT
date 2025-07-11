def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, fungi, and pollen
    to determine the nature of the fungal infections.
    """

    # Baseline data for non-infected bees
    baseline_mortality = 10  # in percent
    
    # --- Pathogenicity Analysis (based on mortality) ---
    print("--- Mortality Analysis ---")
    print(f"A pathogen is defined as an organism causing harm, indicated here by an increased mortality rate.")
    
    # Fungus A data
    fungus_A_mortality = 35 # in percent
    mortality_increase_A = fungus_A_mortality - baseline_mortality
    print(f"Fungus A: Mortality increases by {mortality_increase_A}% (from {baseline_mortality}% to {fungus_A_mortality}%). This indicates it is a pathogen.")

    # Fungus B data
    fungus_B_mortality = 20 # in percent
    mortality_increase_B = fungus_B_mortality - baseline_mortality
    print(f"Fungus B: Mortality increases by {mortality_increase_B}% (from {baseline_mortality}% to {fungus_B_mortality}%). This indicates it is also a pathogen.")

    # Fungus C data
    fungus_C_mortality = 10 # in percent
    mortality_increase_C = fungus_C_mortality - baseline_mortality
    print(f"Fungus C: Mortality increases by {mortality_increase_C}% (from {baseline_mortality}% to {fungus_C_mortality}%). The lack of increased mortality suggests it is not a pathogen.")

    # --- Commensalism Analysis (based on harm/benefit) ---
    print("\n--- Productivity & Commensalism Analysis ---")
    print("A commensal relationship is one where one organism benefits and the other is not significantly harmed.")
    
    # Productivity data for bees fed on Buck pollen
    buck_uninfected_eggs = 45
    buck_infected_A_eggs = 10
    buck_infected_C_eggs = 60
    
    change_A = buck_infected_A_eggs - buck_uninfected_eggs
    print(f"Fungus A (Pathogen): Egg production changes by {change_A} (from {buck_uninfected_eggs} to {buck_infected_A_eggs}), confirming harm to the colony's productivity.")
    
    change_C = buck_infected_C_eggs - buck_uninfected_eggs
    print(f"Fungus C (Commensal): Egg production changes by +{change_C} (from {buck_uninfected_eggs} to {buck_infected_C_eggs}), showing no harm and even a potential benefit.")
    
    print("\nConclusion: The data strongly supports that Fungus A and B are pathogens, while Fungus C is a commensal.")

analyze_bee_experiments()
<<<I>>>