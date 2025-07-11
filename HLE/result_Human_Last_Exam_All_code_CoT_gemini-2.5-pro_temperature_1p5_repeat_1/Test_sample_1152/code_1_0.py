def analyze_bee_experiments():
    """
    Analyzes the bee experiment data to determine the correct statement.
    """

    # --- Data from the experiments ---
    mortality_non_infected = 10  # %
    
    # Experiment 1: Fungus A
    mortality_fungus_A_max = 35 # % (for Buck, Lavender, etc.)
    
    # Experiment 3: Fungus B
    mortality_fungus_B = 20 # % (for all pollens)

    # Experiment 4 & 5: Fungus C
    mortality_fungus_C = 10 # % (for all pollens)
    eggs_buck_uninfected = 45
    eggs_buck_infected_C = 60
    eggs_lavender_uninfected = 30
    eggs_lavender_infected_C = 50

    # --- Step-by-step analysis ---
    
    print("Step 1: Analyze if Fungus A and Fungus B are pathogens.")
    print(f"The baseline mortality rate for non-infected bees is {mortality_non_infected}%.")
    
    # Analysis for Fungus A
    print(f"The highest mortality rate for bees infected with Fungus A is {mortality_fungus_A_max}%.")
    print(f"Comparing Fungus A mortality to baseline: {mortality_fungus_A_max}% > {mortality_non_infected}%.")
    print("Since the mortality rate is higher, Fungus A is a pathogen.")

    # Analysis for Fungus B
    print(f"\nThe highest mortality rate for bees infected with Fungus B is {mortality_fungus_B}%.")
    print(f"Comparing Fungus B mortality to baseline: {mortality_fungus_B}% > {mortality_non_infected}%.")
    print("Since the mortality rate is higher, Fungus B is also a pathogen.")
    
    print("\nConclusion for Step 1: The statement 'Fungus A and B are pathogens' is TRUE.")

    print("\n" + "="*50 + "\n")

    print("Step 2: Analyze the nature of Fungus C.")
    print(f"The mortality rate for bees infected with Fungus C is {mortality_fungus_C}%.")
    print(f"Comparing Fungus C mortality to baseline: {mortality_fungus_C}% is equal to {mortality_non_infected}%.")
    print("Fungus C does not increase the mortality rate, so it does not harm the bees' survival.")

    print("\nNow, let's check the effect of Fungus C on productivity (egg laying):")
    print(f"For bees fed on buck pollen, egg count changed from {eggs_buck_uninfected} (uninfected) to {eggs_buck_infected_C} (infected).")
    print(f"For bees fed on lavender pollen, egg count changed from {eggs_lavender_uninfected} (uninfected) to {eggs_lavender_infected_C} (infected).")
    print("Productivity generally increased or was not negatively impacted in the same way as with pathogens A and B.")
    print("Since Fungus C does not cause harm (no increased mortality) and may provide a benefit (increased productivity), it can be classified as a commensal (or even mutualistic) organism, not a pathogen.")

    print("\nConclusion for Step 2: The statement 'Fungus C is a commensal' is TRUE.")
    
    print("\n" + "="*50 + "\n")
    
    print("Final Conclusion:")
    print("Both parts of statement I are supported by the data:")
    print("1. Fungus A and B are pathogens. (TRUE)")
    print("2. Fungus C is a commensal. (TRUE)")
    print("Therefore, statement I is the correct and most complete answer.")
    
    print("<<<I>>>")

analyze_bee_experiments()