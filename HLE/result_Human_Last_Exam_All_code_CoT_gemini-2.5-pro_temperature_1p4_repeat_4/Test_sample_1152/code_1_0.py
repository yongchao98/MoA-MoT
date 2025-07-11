def analyze_fungi_data():
    """
    Analyzes the provided experimental data to determine the nature of
    Fungus A, B, and C and selects the correct answer choice.
    """

    # Data from the experiments
    mortality_non_infected = 10  # %
    
    # Fungus A data
    mortality_fungus_A = 35  # % (for most pollens)
    
    # Fungus B data
    mortality_fungus_B = 20  # %
    
    # Fungus C data
    mortality_fungus_C = 10  # %
    productivity_C_infected_buck = 60 # eggs
    productivity_C_uninfected_buck = 45 # eggs

    print("Step 1: Analyze if Fungus A and B are pathogens.")
    print("A pathogen causes harm, such as increased mortality.")
    print(f"The baseline mortality rate for non-infected honeybees is {mortality_non_infected}%.")
    
    print("\nAnalyzing Fungus A...")
    print(f"The mortality rate for bees infected with Fungus A is {mortality_fungus_A}%.")
    is_pathogen_A = mortality_fungus_A > mortality_non_infected
    print(f"Is {mortality_fungus_A}% > {mortality_non_infected}%? {is_pathogen_A}.")
    if is_pathogen_A:
        print("Conclusion: Fungus A increases mortality, so it is a pathogen.")
    
    print("\nAnalyzing Fungus B...")
    print(f"The mortality rate for bees infected with Fungus B is {mortality_fungus_B}%.")
    is_pathogen_B = mortality_fungus_B > mortality_non_infected
    print(f"Is {mortality_fungus_B}% > {mortality_non_infected}%? {is_pathogen_B}.")
    if is_pathogen_B:
        print("Conclusion: Fungus B increases mortality, so it is a pathogen.")

    print("\n----------------------------------------------------")
    print("\nStep 2: Analyze if Fungus C is a commensal.")
    print("A commensal organism does not harm its host.")
    print(f"The mortality rate for bees infected with Fungus C is {mortality_fungus_C}%.")
    is_commensal_C = mortality_fungus_C <= mortality_non_infected
    print(f"Is {mortality_fungus_C}% <= {mortality_non_infected}%? {is_commensal_C}.")
    if is_commensal_C:
        print("Conclusion: Fungus C does not increase the mortality rate.")
        print("Furthermore, productivity for bees infected with Fungus C (e.g., producing",
              f"{productivity_C_infected_buck} eggs on buck pollen) can be higher than uninfected bees",
              f"(who produce {productivity_C_uninfected_buck} eggs).")
        print("Since it causes no harm and may even be beneficial, Fungus C is a commensal (or mutualist).")

    print("\n----------------------------------------------------")
    print("\nStep 3: Combine findings and choose the correct answer.")
    print("Based on the analysis:")
    print("1. Fungus A is a pathogen.")
    print("2. Fungus B is a pathogen.")
    print("3. Fungus C is a commensal.")
    print("\nThe statement 'Fungus A and B are pathogens. Fungus C is a commensal.' accurately summarizes the data.")
    print("This corresponds to answer choice I.")

# Execute the analysis
analyze_fungi_data()
