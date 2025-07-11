import sys
# Helper function to print to stderr for thinking process
def thinking_print(s):
    sys.stderr.write(s + "\n")

# --- Data from the experiments ---
baseline_mortality = 0.10 # 10%

# Experiment 1 & 2: Fungus A
fungus_a_data = {
    'mortality_rates': [0.35, 0.10, 0.35, 0.35, 0.35, 0.35, 0.35] # buck, sunflower, lavender, etc.
}

# Experiment 3: Fungus B
fungus_b_data = {
    'mortality_rates': [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20]
}

# Experiment 4 & 5: Fungus C
fungus_c_data = {
    'mortality_rates': [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
    'productivity_infected': {
        "buck": 60, "sunflower": 25, "lavender": 50, "canola": 50,
        "milkweed": 50, "aster": 50, "mixed": 52
    },
    'productivity_uninfected': {
        "buck": 45, "sunflower": 30, "lavender": 30, "canola": 30,
        "milkweed": 30, "aster": 30, "mixed": 32
    }
}


def analyze_fungi():
    """
    Analyzes the experimental data to determine the nature of each fungus.
    """
    print("Step 1: Analyzing Fungus A")
    max_mortality_A = max(fungus_a_data['mortality_rates'])
    print(f"The maximum mortality rate for bees infected with Fungus A is {max_mortality_A * 100}%.")
    print(f"The baseline mortality rate for non-infected bees is {baseline_mortality * 100}%.")
    if max_mortality_A > baseline_mortality:
        increase = (max_mortality_A - baseline_mortality) * 100
        print(f"Since {max_mortality_A * 100}% > {baseline_mortality * 100}%, Fungus A increases mortality. It is a pathogen.")
    else:
        print("Fungus A does not increase mortality. It is not a pathogen.")

    print("\nStep 2: Analyzing Fungus B")
    max_mortality_B = max(fungus_b_data['mortality_rates'])
    print(f"The maximum mortality rate for bees infected with Fungus B is {max_mortality_B * 100}%.")
    print(f"The baseline mortality rate is {baseline_mortality * 100}%.")
    if max_mortality_B > baseline_mortality:
        increase = (max_mortality_B - baseline_mortality) * 100
        print(f"Since {max_mortality_B * 100}% > {baseline_mortality * 100}%, Fungus B increases mortality. It is a pathogen.")
    else:
        print("Fungus B does not increase mortality. It is not a pathogen.")

    print("\nStep 3: Analyzing Fungus C")
    max_mortality_C = max(fungus_c_data['mortality_rates'])
    print(f"The maximum mortality rate for bees infected with Fungus C is {max_mortality_C * 100}%.")
    print(f"The baseline mortality rate is {baseline_mortality * 100}%.")
    if max_mortality_C > baseline_mortality:
        print("Fungus C increases mortality. It is a pathogen.")
    else:
        print(f"Since {max_mortality_C * 100}% is not greater than {baseline_mortality * 100}%, Fungus C does not increase mortality.")
        print("Checking productivity effects:")
        # Check if productivity generally increases
        increases = 0
        for pollen in fungus_c_data['productivity_infected']:
            infected_eggs = fungus_c_data['productivity_infected'][pollen]
            uninfected_eggs = fungus_c_data['productivity_uninfected'][pollen]
            if infected_eggs > uninfected_eggs:
                increases += 1
                thinking_print(f"  - For {pollen} pollen, productivity increased from {uninfected_eggs} to {infected_eggs} eggs.")
        print(f"Productivity (egg laying) increased in {increases} out of 7 cases.")
        print("Because Fungus C does not cause harm (increased mortality) and can provide a benefit (increased productivity), it can be classified as a commensal or mutualist.")

    print("\nConclusion:")
    print("Based on the analysis, Fungus A and B are pathogens, and Fungus C is a commensal. This corresponds to answer choice I.")


# Run the analysis
analyze_fungi()
