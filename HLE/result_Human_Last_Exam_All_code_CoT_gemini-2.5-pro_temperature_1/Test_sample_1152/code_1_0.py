def analyze_bee_experiments():
    """
    Analyzes the provided experimental data to determine the nature of each fungus
    and select the correct conclusion from the given options.
    """

    # --- Data from the problem description ---
    control_mortality = 10  # Mortality rate of non-infected honeybees in %

    # Fungus A Data (untreated by special pollen, e.g., from buck)
    fungus_A_mortality = 35  # %

    # Fungus B Data
    fungus_B_mortality = 20  # %

    # Fungus C Data
    fungus_C_mortality = 10  # %
    # Productivity data for Fungus C (using buck pollen as an example)
    fungus_C_eggs_uninfected_buck = 45
    fungus_C_eggs_infected_buck = 60


    # --- Step-by-step Analysis ---

    print("Step 1: Analyze Fungus A")
    print(f"The mortality rate for bees infected with Fungus A is {fungus_A_mortality}%.")
    print(f"The mortality rate for non-infected bees is {control_mortality}%.")
    is_A_pathogen = fungus_A_mortality > control_mortality
    print(f"Conclusion: Since {fungus_A_mortality} > {control_mortality}, Fungus A increases mortality and is a pathogen.\n")


    print("Step 2: Analyze Fungus B")
    print(f"The mortality rate for bees infected with Fungus B is {fungus_B_mortality}%.")
    print(f"The mortality rate for non-infected bees is {control_mortality}%.")
    is_B_pathogen = fungus_B_mortality > control_mortality
    print(f"Conclusion: Since {fungus_B_mortality} > {control_mortality}, Fungus B increases mortality and is a pathogen.\n")


    print("Step 3: Analyze Fungus C")
    print(f"The mortality rate for bees infected with Fungus C is {fungus_C_mortality}%.")
    print(f"The mortality rate for non-infected bees is {control_mortality}%.")
    is_C_pathogen = fungus_C_mortality > control_mortality
    if not is_C_pathogen:
        print(f"Conclusion: Since the mortality rate {fungus_C_mortality}% is equal to the control rate {control_mortality}%, Fungus C does not cause additional death.")
        # Check productivity data
        productivity_increase = fungus_C_eggs_infected_buck > fungus_C_eggs_uninfected_buck
        print(f"Furthermore, with Fungus C, egg production for bees fed buck pollen increased from {fungus_C_eggs_uninfected_buck} to {fungus_C_eggs_infected_buck}.")
        print("An organism that does not cause harm and may provide a benefit is a commensal or mutualist.\n")
    else:
        print("Conclusion: Fungus C is a pathogen.\n")


    # --- Final Conclusion ---
    print("Summary of Findings:")
    print(f"1. Fungus A is a pathogen: {is_A_pathogen}")
    print(f"2. Fungus B is a pathogen: {is_B_pathogen}")
    print(f"3. Fungus C is a commensal: {not is_C_pathogen}")
    print("\nBased on this analysis, the statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is the most accurate conclusion.")

if __name__ == "__main__":
    analyze_bee_experiments()
<<<I>>>