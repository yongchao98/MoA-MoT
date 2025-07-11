import pandas as pd

def analyze_bee_data():
    """
    Analyzes honeybee experimental data to determine the correct conclusion.
    """
    # Baseline data
    non_infected_mortality = 10  # %

    # --- Data Organization ---
    # Fungus A Data (Exp 1 & 2)
    fungus_A_data = {
        'Pollen': ['Buck', 'Sunflower', 'Lavender', 'Canola', 'Milkweed', 'Aster', 'Mixed'],
        'Fungal Cells (per 0.02ul)': [200, 35, 200, 200, 200, 200, 180],
        'Mortality (%)': [35, 10, 35, 35, 35, 35, 35],
        'Eggs (Not Infected)': [45, 30, 30, 30, 30, 30, 32],
        'Eggs (Infected)': [10, 20, 10, 8, 9, 11, 12]
    }

    # Fungus B Data (Exp 3)
    fungus_B_data = {
        'Pollen': ['Buck', 'Sunflower', 'Lavender', 'Canola', 'Milkweed', 'Aster', 'Mixed'],
        'Mortality (%)': [20, 20, 20, 20, 20, 20, 20]
    }

    # Fungus C Data (Exp 4 & 5)
    fungus_C_data = {
        'Pollen': ['Buck', 'Sunflower', 'Lavender', 'Canola', 'Milkweed', 'Aster', 'Mixed'],
        'Mortality (%)': [10, 10, 10, 10, 10, 10, 10],
        'Eggs (Infected)': [60, 25, 50, 50, 50, 50, 52]
    }

    print("Step 1: Determine which fungi are pathogens by checking if they increase mortality.")
    # Fungus A mortality (most pollens)
    mortality_A = fungus_A_data['Mortality (%)'][0]
    print(f"Non-infected mortality is {non_infected_mortality}%.")
    print(f"Fungus A mortality is {mortality_A}%. Since {mortality_A} > {non_infected_mortality}, Fungus A is a pathogen.")

    # Fungus B mortality
    mortality_B = fungus_B_data['Mortality (%)'][0]
    print(f"Fungus B mortality is {mortality_B}%. Since {mortality_B} > {non_infected_mortality}, Fungus B is a pathogen.")
    
    # Fungus C mortality
    mortality_C = fungus_C_data['Mortality (%)'][0]
    print(f"Fungus C mortality is {mortality_C}%. Since {mortality_C} == {non_infected_mortality}, Fungus C is not a pathogen.\n")

    print("Step 2: Compare the deadliness and treatability of Fungus A and Fungus B.")
    print(f"Comparing mortality rates: Fungus A ({mortality_A}%) vs Fungus B ({mortality_B}%).")
    print(f"Since {mortality_A} > {mortality_B}, Fungus A is more deadly than Fungus B.")

    mortality_A_sunflower = fungus_A_data['Mortality (%)'][1]
    mortality_B_sunflower = fungus_B_data['Mortality (%)'][1]
    print(f"With sunflower pollen, Fungus A mortality drops from {mortality_A}% to {mortality_A_sunflower}%.")
    print(f"With sunflower pollen, Fungus B mortality remains at {mortality_B_sunflower}%.")
    print("Therefore, Fungus A is easier to treat (with sunflower pollen) than Fungus B, making Fungus B more difficult to treat.\n")

    print("Step 3: Evaluate claims about buck pollen and productivity.")
    eggs_buck_uninfected = fungus_A_data['Eggs (Not Infected)'][0]
    eggs_sunflower_uninfected = fungus_A_data['Eggs (Not Infected)'][1]
    print(f"When not infected, bees fed buck pollen produce {eggs_buck_uninfected} eggs, while those on sunflower produce {eggs_sunflower_uninfected}.")

    eggs_buck_infected_A = fungus_A_data['Eggs (Infected)'][0]
    eggs_sunflower_infected_A = fungus_A_data['Eggs (Infected)'][1]
    print(f"When infected with Fungus A, bees fed buck pollen produce {eggs_buck_infected_A} eggs, while those on sunflower produce {eggs_sunflower_infected_A}.")
    print(f"Since {eggs_buck_infected_A} < {eggs_sunflower_infected_A}, the claim that buck pollen *always* leads to more eggs is false.")
    print("The statement 'Honeybees fed on buck will not always produce more eggs than honeybees fed on any other pollen' is true.\n")

    print("Step 4: Final Conclusion based on evaluating the answer choices.")
    print("Let's evaluate option F: 'Fungus A is more deadly than fungus B. Fungus B is more difficult to treat than fungus A. Honeybees fed on buck will not always produce more eggs than honeybees fed on any other pollen.'")
    print(f" - 'Fungus A is more deadly than fungus B' is TRUE ({mortality_A}% > {mortality_B}%).")
    print(f" - 'Fungus B is more difficult to treat than fungus A' is TRUE (sunflower pollen only reduces mortality for A).")
    print(f" - 'Honeybees fed on buck will not always produce more eggs...' is TRUE ({eggs_buck_infected_A} eggs < {eggs_sunflower_infected_A} eggs).")
    print("\nAll parts of statement F are correct and supported by the data.")

analyze_bee_data()
<<<F>>>