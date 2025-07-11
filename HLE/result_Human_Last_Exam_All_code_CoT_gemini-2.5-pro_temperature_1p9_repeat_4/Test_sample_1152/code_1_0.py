import collections

def analyze_bee_experiments():
    """
    Analyzes data from five bee experiments to draw conclusions.
    """
    # --- Data from Experiments ---

    # Baseline for healthy bees
    non_infected_mortality = 10  # percent

    # Experiment 1 & 2: Fungus A (Intestine)
    exp1_fungus_a_mortality = {
        'Buck': 35, 'Sunflower': 10, 'Lavender': 35, 'Canola': 35,
        'Milkweed': 35, 'Aster': 35, 'Mixed': 35
    }
    exp2_fungus_a_productivity = {
        'Buck': {'not_infected': 45, 'infected': 10},
        'Sunflower': {'not_infected': 30, 'infected': 20},
        'Lavender': {'not_infected': 30, 'infected': 10},
        'Canola': {'not_infected': 30, 'infected': 8},
        'Milkweed': {'not_infected': 30, 'infected': 9},
        'Aster': {'not_infected': 30, 'infected': 11},
        'Mixed': {'not_infected': 32, 'infected': 12}
    }

    # Experiment 3: Fungus B (Surface)
    exp3_fungus_b_mortality = {
        'Buck': 20, 'Sunflower': 20, 'Lavender': 20, 'Canola': 20,
        'Milkweed': 20, 'Aster': 20, 'Mixed': 20
    }

    # Experiment 4 & 5: Fungus C (Intestine)
    exp4_fungus_c_mortality = {
        'Buck': 10, 'Sunflower': 10, 'Lavender': 10, 'Canola': 10,
        'Milkweed': 10, 'Aster': 10, 'Mixed': 10
    }
    exp5_fungus_c_productivity = {
        'Buck': {'not_infected': 45, 'infected': 60},
        'Sunflower': {'not_infected': 30, 'infected': 25},
        'Lavender': {'not_infected': 30, 'infected': 50},
        'Canola': {'not_infected': 30, 'infected': 50},
        'Milkweed': {'not_infected': 30, 'infected': 50},
        'Aster': {'not_infected': 30, 'infected': 50},
        'Mixed': {'not_infected': 32, 'infected': 52}
    }

    print("--- Analysis Step 1: Classify Fungi based on Mortality ---")

    # Fungus A Analysis
    max_mortality_A = max(exp1_fungus_a_mortality.values())
    if max_mortality_A > non_infected_mortality:
        print(f"Fungus A is a pathogen: Its max mortality rate ({max_mortality_A}%) is higher than the non-infected rate ({non_infected_mortality}%).")
    else:
        print("Fungus A is not a pathogen.")

    # Fungus B Analysis
    max_mortality_B = max(exp3_fungus_b_mortality.values())
    if max_mortality_B > non_infected_mortality:
        print(f"Fungus B is a pathogen: Its max mortality rate ({max_mortality_B}%) is higher than the non-infected rate ({non_infected_mortality}%).")
    else:
        print("Fungus B is not a pathogen.")

    # Fungus C Analysis
    max_mortality_C = max(exp4_fungus_c_mortality.values())
    if max_mortality_C > non_infected_mortality:
        print(f"Fungus C is a pathogen.")
    else:
        # Check for productivity increase with Fungus C
        prod_increase = exp5_fungus_c_productivity['Buck']['infected'] > exp5_fungus_c_productivity['Buck']['not_infected']
        print(f"Fungus C is a commensal/mutualist: Its mortality rate ({max_mortality_C}%) is equal to the non-infected rate ({non_infected_mortality}%).")
        print(f"It can even be beneficial, increasing productivity (e.g., Buck pollen: {exp5_fungus_c_productivity['Buck']['infected']} eggs infected vs. {exp5_fungus_c_productivity['Buck']['not_infected']} eggs not infected).")

    print("\n--- Analysis Step 2: Compare Pathogens A and B ---")
    print(f"Deadliness: Fungus A (max mortality {max_mortality_A}%) is more deadly than Fungus B (max mortality {max_mortality_B}%).")

    sunflower_mortality_A = exp1_fungus_a_mortality['Sunflower']
    print(f"Treatability of A: Fungus A is treatable. Sunflower pollen reduces mortality from {max_mortality_A}% to {sunflower_mortality_A}%.")
    print(f"Treatability of B: Fungus B appears untreatable with tested pollens, as mortality remains at {max_mortality_B}%.")
    print("Conclusion: Fungus B is more difficult to treat than Fungus A.")


    print("\n--- Analysis Step 3: Evaluate Pollen Claims ---")
    # Check "Buck always produces more eggs"
    buck_infected_A = exp2_fungus_a_productivity['Buck']['infected']
    sunflower_infected_A = exp2_fungus_a_productivity['Sunflower']['infected']
    print("Is Buck pollen always best for productivity? No.")
    print(f"Counterexample: For bees with Fungus A, Sunflower is more productive ({sunflower_infected_A} eggs) than Buck ({buck_infected_A} eggs).")


    print("\n--- FINAL CONCLUSION ---")
    print("The data shows that Fungus A and B are pathogens because they increase mortality above the baseline.")
    print("Fungus C is commensal because it does not increase mortality and can increase productivity.")
    print("This directly supports the statement: 'Fungus A and B are pathogens. Fungus C is a commensal.'")

analyze_bee_experiments()