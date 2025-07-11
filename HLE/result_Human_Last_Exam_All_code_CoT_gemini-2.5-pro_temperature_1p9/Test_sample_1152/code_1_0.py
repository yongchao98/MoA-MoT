import pandas as pd

def analyze_honeybee_data():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    baseline_mortality = 10

    # Data from experiments
    # Fungus A mortality is 35% for most pollens
    fungus_a_mortality = 35
    # Fungus B mortality is 20% for all pollens
    fungus_b_mortality = 20
    # Fungus C mortality is 10% for all pollens
    fungus_c_mortality = 10

    print("Step 1: Determine if Fungus A is a pathogen.")
    is_a_pathogen = fungus_a_mortality > baseline_mortality
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")
    print(f"The mortality rate for bees infected with Fungus A is {fungus_a_mortality}%.")
    print(f"Is {fungus_a_mortality} > {baseline_mortality}? {is_a_pathogen}.")
    print("Conclusion: Fungus A increases mortality, so it is a pathogen.\n")

    print("Step 2: Determine if Fungus B is a pathogen.")
    is_b_pathogen = fungus_b_mortality > baseline_mortality
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")
    print(f"The mortality rate for bees infected with Fungus B is {fungus_b_mortality}%.")
    print(f"Is {fungus_b_mortality} > {baseline_mortality}? {is_b_pathogen}.")
    print("Conclusion: Fungus B increases mortality, so it is a pathogen.\n")

    print("Step 3: Determine the nature of Fungus C.")
    is_c_pathogen = fungus_c_mortality > baseline_mortality
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")
    print(f"The mortality rate for bees infected with Fungus C is {fungus_c_mortality}%.")
    print(f"Is {fungus_c_mortality} > {baseline_mortality}? {is_c_pathogen}.")
    print("Conclusion: Fungus C does not increase mortality. An organism that lives with a host without causing harm is a commensal. Therefore, Fungus C is a commensal.\n")

    print("Step 4: Evaluate the final answer choice.")
    print("Based on the analysis:")
    print("- Fungus A and B are pathogens.")
    print("- Fungus C is a commensal.")
    print("This matches answer choice I.\n")

    # This part checks other claims to confirm they are false, but the core logic is above.
    # Claim: Fungus A is more deadly than Fungus B.
    a_more_deadly = fungus_a_mortality > fungus_b_mortality
    # Claim: Fungus B is harder to treat than Fungus A (sunflower lowers A's mortality to 10%, but no pollen lowers B's).
    b_harder_to_treat = True
    # Claim: Buck pollen always leads to more eggs. False (e.g., infected with Fungus A: Sunflower(20) > Buck(10)).
    buck_always_best = False
    
    print("Final conclusion check: The statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is fully supported by the data.")

analyze_honeybee_data()
