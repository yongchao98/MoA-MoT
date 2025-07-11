import pandas as pd

def analyze_honeybee_data():
    """
    Analyzes the provided honeybee experimental data to draw conclusions
    and select the correct statement.
    """
    # --- Data Representation ---
    control_mortality = 10  # %

    # Experiment 1 & 2: Fungus A Data
    fungus_A = {
        'buck':      {'fungal_cells_0.02ul': 200, 'mortality': 35, 'eggs_not_infected': 45, 'eggs_infected': 10},
        'sunflower': {'fungal_cells_0.02ul': 35,  'mortality': 10, 'eggs_not_infected': 30, 'eggs_infected': 20},
        'lavender':  {'fungal_cells_0.02ul': 200, 'mortality': 35, 'eggs_not_infected': 30, 'eggs_infected': 10},
        'canola':    {'fungal_cells_0.02ul': 200, 'mortality': 35, 'eggs_not_infected': 30, 'eggs_infected': 8},
        'milkweed':  {'fungal_cells_0.02ul': 200, 'mortality': 35, 'eggs_not_infected': 30, 'eggs_infected': 9},
        'aster':     {'fungal_cells_0.02ul': 200, 'mortality': 35, 'eggs_not_infected': 30, 'eggs_infected': 11}
    }

    # Experiment 3: Fungus B Data
    fungus_B = {
        'buck':      {'mortality': 20},
        'sunflower': {'mortality': 20},
        # ... other pollens have the same 20% mortality
    }
    
    # Experiment 4 & 5: Fungus C Data
    fungus_C = {
        'buck':      {'mortality': 10, 'eggs_not_infected': 45, 'eggs_infected': 60},
        'sunflower': {'mortality': 10, 'eggs_not_infected': 30, 'eggs_infected': 25},
        'lavender':  {'mortality': 10, 'eggs_not_infected': 30, 'eggs_infected': 50},
    }

    print("--- Analysis of Experimental Data ---")

    # 1. Analyze Pathogenicity of Fungus A and B
    print("\nStep 1: Analyzing Pathogenicity")
    mortality_A_untreated = fungus_A['buck']['mortality']
    print(f"Is Fungus A a pathogen? Comparing its mortality rate to control:")
    print(f"  - Mortality with Fungus A (e.g., on buck pollen): {mortality_A_untreated}%")
    print(f"  - Control mortality (not infected): {control_mortality}%")
    print(f"  - Conclusion: Since {mortality_A_untreated} > {control_mortality}, Fungus A is a pathogen.")

    mortality_B = fungus_B['buck']['mortality']
    print(f"\nIs Fungus B a pathogen? Comparing its mortality rate to control:")
    print(f"  - Mortality with Fungus B: {mortality_B}%")
    print(f"  - Control mortality (not infected): {control_mortality}%")
    print(f"  - Conclusion: Since {mortality_B} > {control_mortality}, Fungus B is a pathogen.")

    # 2. Analyze Pathogenicity and Nature of Fungus C
    print("\nStep 2: Analyzing Fungus C")
    mortality_C = fungus_C['buck']['mortality']
    print(f"Is Fungus C a pathogen? Comparing its mortality rate to control:")
    print(f"  - Mortality with Fungus C: {mortality_C}%")
    print(f"  - Control mortality (not infected): {control_mortality}%")
    print(f"  - Conclusion: Since {mortality_C} == {control_mortality}, Fungus C does not increase mortality.")
    
    eggs_C_not_infected = fungus_C['buck']['eggs_not_infected']
    eggs_C_infected = fungus_C['buck']['eggs_infected']
    print(f"\nDoes Fungus C harm productivity?")
    print(f"  - Egg production with Fungus C (buck pollen): {eggs_C_infected}")
    print(f"  - Egg production without infection (buck pollen): {eggs_C_not_infected}")
    print(f"  - Conclusion: Since {eggs_C_infected} > {eggs_C_not_infected}, Fungus C appears to increase productivity, suggesting it is not a pathogen. The term 'commensal' (or mutualist) is appropriate.")
    
    print("\n--- Final Summary based on Analysis ---")
    print("1. Fungus A is a pathogen (increases mortality).")
    print("2. Fungus B is a pathogen (increases mortality).")
    print("3. Fungus C is not a pathogen; it doesn't increase mortality and can increase productivity, so it can be described as a commensal.")

analyze_honeybee_data()
print("\nBased on the analysis, the most accurate and comprehensive statement is that Fungus A and B are pathogens, and Fungus C is a commensal.")
<<<I>>>