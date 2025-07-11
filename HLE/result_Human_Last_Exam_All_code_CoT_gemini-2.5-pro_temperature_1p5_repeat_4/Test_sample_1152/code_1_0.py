def analyze_bee_data():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    # Baseline data
    control_mortality = 10

    # Data from experiments for each fungus
    # Fungus A (pathogen)
    fungus_A_max_mortality = 35
    # Fungus B (pathogen)
    fungus_B_max_mortality = 20
    # Fungus C (commensal)
    fungus_C_max_mortality = 10
    # Productivity data for Fungus C
    eggs_not_infected_buck = 45
    eggs_infected_C_buck = 60

    print("Step 1: Analyze Fungus A's effect on mortality.")
    print(f"The baseline mortality for non-infected bees is {control_mortality}%.")
    print(f"The maximum mortality for bees infected with Fungus A is {fungus_A_max_mortality}%.")
    is_pathogen_A = fungus_A_max_mortality > control_mortality
    print(f"Is Fungus A a pathogen? ({fungus_A_max_mortality} > {control_mortality}) -> {is_pathogen_A}. Fungus A is a pathogen.")
    print("-" * 20)

    print("Step 2: Analyze Fungus B's effect on mortality.")
    print(f"The baseline mortality for non-infected bees is {control_mortality}%.")
    print(f"The mortality for bees infected with Fungus B is {fungus_B_max_mortality}%.")
    is_pathogen_B = fungus_B_max_mortality > control_mortality
    print(f"Is Fungus B a pathogen? ({fungus_B_max_mortality} > {control_mortality}) -> {is_pathogen_B}. Fungus B is a pathogen.")
    print("-" * 20)
    
    print("Step 3: Analyze Fungus C's effect on mortality and productivity.")
    print(f"The baseline mortality for non-infected bees is {control_mortality}%.")
    print(f"The mortality for bees infected with Fungus C is {fungus_C_max_mortality}%.")
    is_pathogen_C = fungus_C_max_mortality > control_mortality
    print(f"Is Fungus C a pathogen? ({fungus_C_max_mortality} > {control_mortality}) -> {is_pathogen_C}. Fungus C does not increase mortality.")
    
    print("\nAdditionally, let's check its effect on productivity (egg-laying).")
    print(f"For bees on buck pollen, non-infected colonies produced {eggs_not_infected_buck} eggs, while colonies infected with Fungus C produced {eggs_infected_C_buck} eggs.")
    print("Since infection did not increase mortality and in many cases increased productivity, Fungus C can be classified as a commensal (or even mutualist).")
    print("-" * 20)

    print("\nFinal Conclusion:")
    print("The data shows that Fungus A and Fungus B are pathogens because they increase the mortality rate above the baseline.")
    print("Fungus C is a commensal because it does not increase the mortality rate and even shows benefits to productivity.")
    print("Therefore, the correct statement is: Fungus A and B are pathogens. Fungus C is a commensal.")

analyze_bee_data()
<<<I>>>