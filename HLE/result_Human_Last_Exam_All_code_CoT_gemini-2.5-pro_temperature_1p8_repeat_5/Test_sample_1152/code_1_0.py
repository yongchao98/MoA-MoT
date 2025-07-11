import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string
sys.stdout = io.StringIO()

def analyze_honeybee_data():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct conclusion.
    """
    # Baseline data
    baseline_mortality = 10  # % mortality for non-infected bees

    # Fungi mortality data from experiments 1, 3, 4
    # For each fungus, we take the highest observed mortality rate.
    mortality_A = 35  # Highest mortality rate with Fungus A infection
    mortality_B = 20  # Highest mortality rate with Fungus B infection
    mortality_C = 10  # Highest mortality rate with Fungus C infection

    # Productivity data for Fungus C from experiment 5
    # Comparing not infected vs. infected with Fungus C for Buck pollen
    productivity_buck_uninfected = 45
    productivity_buck_infected_C = 60
    
    # --- Analysis ---
    
    # Step 1: Analyze Fungus A
    print("--- Analysis of Fungus A ---")
    print(f"A pathogen is defined as an organism that increases host mortality.")
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")
    print(f"The mortality rate for bees infected with Fungus A is up to {mortality_A}%.")
    print(f"Comparison: {mortality_A}% > {baseline_mortality}%")
    print("Conclusion: Fungus A is a pathogen.\n")

    # Step 2: Analyze Fungus B
    print("--- Analysis of Fungus B ---")
    print(f"The mortality rate for bees infected with Fungus B is {mortality_B}%.")
    print(f"Comparison: {mortality_B}% > {baseline_mortality}%")
    print("Conclusion: Fungus B is a pathogen.\n")
    
    # Step 3: Analyze Fungus C
    print("--- Analysis of Fungus C ---")
    print(f"A commensal organism does not harm the host.")
    print(f"The mortality rate for bees infected with Fungus C is {mortality_C}%.")
    print(f"Comparison: {mortality_C}% is not greater than the baseline {baseline_mortality}% mortality.")
    print(f"This indicates Fungus C is not pathogenic.")
    
    print(f"Furthermore, productivity with Fungus C infection increased from {productivity_buck_uninfected} to {productivity_buck_infected_C} eggs in one case.")
    print("Conclusion: Fungus C is a commensal (or potentially mutualistic).\n")

    # Final Conclusion based on analysis
    print("--- Final Verdict ---")
    print("The evidence shows that Fungus A and B are pathogens because they increase mortality.")
    print("Fungus C is a commensal because it does not increase mortality and may even increase productivity.")
    print("This directly supports the statement: 'Fungus A and B are pathogens. Fungus C is a commensal.'")

# Run analysis
analyze_honeybee_data()

# Get the output and restore stdout
output = sys.stdout.getvalue()
sys.stdout = stdout_backup

# Print the captured output
print(output)
print("<<<I>>>")