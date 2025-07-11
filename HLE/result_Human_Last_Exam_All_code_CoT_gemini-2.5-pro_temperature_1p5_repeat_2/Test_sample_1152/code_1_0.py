import pandas as pd

def analyze_bee_data():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the most accurate conclusion.
    """

    # Data points
    baseline_mortality = 10

    # Fungus A data
    fungus_a_mortality_general = 35
    fungus_a_mortality_sunflower = 10

    # Fungus B data
    fungus_b_mortality = 20

    # Fungus C data
    fungus_c_mortality_general = 10
    fungus_c_productivity_infected_buck = 60
    fungus_c_productivity_not_infected_buck = 45

    # --- Analysis Step 1: Identify Pathogens ---
    print("--- Pathogenicity Analysis (based on mortality rate) ---")
    print(f"Baseline mortality rate of non-infected bees: {baseline_mortality}%")
    print(f"Fungus A causes a mortality rate of up to {fungus_a_mortality_general}%, which is higher than the baseline. Therefore, Fungus A is a pathogen.")
    print(f"Fungus B causes a mortality rate of {fungus_b_mortality}%, which is higher than the baseline. Therefore, Fungus B is a pathogen.")
    
    # --- Analysis Step 2: Characterize Fungus C ---
    print("\n--- Characterization of Fungus C ---")
    print(f"Fungus C causes a mortality rate of {fungus_c_mortality_general}%, which is equal to the non-infected baseline. It does not increase mortality.")
    print(f"Furthermore, Fungus C can increase productivity. E.g., for Buck pollen:")
    print(f" - Eggs in infected colony: {fungus_c_productivity_infected_buck}")
    print(f" - Eggs in not infected colony: {fungus_c_productivity_not_infected_buck}")
    print("Since Fungus C does not cause harm and may be beneficial, it can be classified as a commensal or mutualistic organism.")

    # --- Final Conclusion ---
    print("\n--- Final Conclusion ---")
    print("Based on the analysis, the most accurate statement is:")
    print("'Fungus A and B are pathogens. Fungus C is a commensal.'")


analyze_bee_data()