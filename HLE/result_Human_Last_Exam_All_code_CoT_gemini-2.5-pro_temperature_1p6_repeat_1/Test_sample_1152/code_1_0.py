import pandas as pd

def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to evaluate several hypotheses.
    """

    # --- Data Setup ---
    baseline_mortality = 10  # % for non-infected bees

    data = {
        'Fungus A': {'mortality': 35, 'treatment_mortality': 10},
        'Fungus B': {'mortality': 20, 'treatment_mortality': 20},
        'Fungus C': {'mortality': 10, 'treatment_mortality': 10},
    }

    productivity_data = {
        'Fungus A': {'Buck': (45, 10), 'Sunflower': (30, 20)},
        'Fungus C': {'Buck': (45, 60), 'Sunflower': (30, 25)}
    }

    print("Step-by-step analysis of key claims based on the experimental data:\n")

    # --- Claim 1 & 2: Are Fungus A and B pathogens? ---
    print("1. Determining if Fungus A and B are pathogens...")
    fungus_a_untreated_mortality = data['Fungus A']['mortality']
    is_a_pathogen = fungus_a_untreated_mortality > baseline_mortality
    print(f"   - Fungus A mortality rate ({fungus_a_untreated_mortality}%) vs. Baseline ({baseline_mortality}%).")
    print(f"   - Equation: {fungus_a_untreated_mortality} > {baseline_mortality} is {is_a_pathogen}.")
    print(f"   - Conclusion: Fungus A is a pathogen.\n")

    fungus_b_untreated_mortality = data['Fungus B']['mortality']
    is_b_pathogen = fungus_b_untreated_mortality > baseline_mortality
    print(f"   - Fungus B mortality rate ({fungus_b_untreated_mortality}%) vs. Baseline ({baseline_mortality}%).")
    print(f"   - Equation: {fungus_b_untreated_mortality} > {baseline_mortality} is {is_b_pathogen}.")
    print(f"   - Conclusion: Fungus B is a pathogen.\n")


    # --- Claim 3: Is Fungus C a commensal? ---
    print("2. Determining the nature of Fungus C...")
    fungus_c_mortality = data['Fungus C']['mortality']
    is_c_pathogen = fungus_c_mortality > baseline_mortality
    print(f"   - Fungus C mortality rate ({fungus_c_mortality}%) vs. Baseline ({baseline_mortality}%).")
    print(f"   - Equation: {fungus_c_mortality} > {baseline_mortality} is {is_c_pathogen}.")
    print(f"   - Conclusion: Fungus C does not increase mortality.\n")
    
    # Check productivity impact of Fungus C
    c_buck_healthy, c_buck_infected = productivity_data['Fungus C']['Buck']
    c_sunflower_healthy, c_sunflower_infected = productivity_data['Fungus C']['Sunflower']
    print(f"   - Effect of Fungus C on productivity:")
    print(f"     - With Buck pollen, egg count changes from {c_buck_healthy} to {c_buck_infected}.")
    print(f"     - With Sunflower pollen, egg count changes from {c_sunflower_healthy} to {c_sunflower_infected}.")
    print(f"   - Conclusion: Fungus C does not harm the bees and can even increase productivity, making 'commensal' a fitting description.\n")
    
    print("-" * 50)
    print("Summary of Facts:")
    print("1. Fungus A is a pathogen.")
    print("2. Fungus B is a pathogen.")
    print("3. Fungus C is not a pathogen (commensal or mutualist).")
    print("\nBased on this analysis, the statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is correct.")


analyze_bee_experiments()