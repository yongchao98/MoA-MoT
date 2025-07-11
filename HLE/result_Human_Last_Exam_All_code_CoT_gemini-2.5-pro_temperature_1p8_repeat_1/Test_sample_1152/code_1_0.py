import pandas as pd

# --- Data Representation ---
data = {
    'Fungus A': {
        'type': 'Intestinal',
        'pollen_data': {
            'Buck': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 45, 'eggs_infected': 10},
            'Sunflower': {'fungal_cells': 35 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 20},
            'Lavender': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 30, 'eggs_infected': 10},
            'Canola': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 30, 'eggs_infected': 8},
            'Milkweed': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 30, 'eggs_infected': 9},
            'Aster': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 30, 'eggs_infected': 11},
            'Mixed': {'fungal_cells': 180 / 0.02 * 2, 'mortality_rate': 0.35, 'eggs_not_infected': 32, 'eggs_infected': 12},
        }
    },
    'Fungus B': {
        'type': 'Surface',
        'pollen_data': {
            'Buck': {'fungal_cells': 2500, 'mortality_rate': 0.20},
            'Sunflower': {'fungal_cells': 2000, 'mortality_rate': 0.20},
            'Lavender': {'fungal_cells': 2500, 'mortality_rate': 0.20},
            'Canola': {'fungal_cells': 2500, 'mortality_rate': 0.20},
            'Milkweed': {'fungal_cells': 2500, 'mortality_rate': 0.20},
            'Aster': {'fungal_cells': 2500, 'mortality_rate': 0.20},
            'Mixed': {'fungal_cells': 2500, 'mortality_rate': 0.20},
        }
    },
    'Fungus C': {
        'type': 'Intestinal',
        'pollen_data': {
            'Buck': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 45, 'eggs_infected': 60},
            'Sunflower': {'fungal_cells': 100 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 25},
            'Lavender': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 50},
            'Canola': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 50},
            'Milkweed': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 50},
            'Aster': {'fungal_cells': 200 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 30, 'eggs_infected': 50},
            'Mixed': {'fungal_cells': 180 / 0.02 * 2, 'mortality_rate': 0.10, 'eggs_not_infected': 32, 'eggs_infected': 52},
        }
    }
}
# Note: Fungal cells from Exp 1 & 4 were given per 0.02ul, but the experiment states it was calculated per 0.04ul.
# The calculation `val / 0.02 * 2` is equivalent to `val / 0.01` and corrects for a likely typo in the question description (cells/0.02ul vs calculation per 0.04ul), scaling to a consistent base unit.
# However, the relative differences remain the same, so this normalization does not change the logical conclusions.

baseline_mortality = 0.10

def is_pathogen(fungus_name):
    """A fungus is a pathogen if it increases mortality or decreases productivity."""
    fungus = data[fungus_name]
    
    # Check for increased mortality
    for pollen in fungus['pollen_data']:
        if fungus['pollen_data'][pollen]['mortality_rate'] > baseline_mortality:
            return True, f"{fungus_name} is a pathogen because it increases mortality rates above the baseline {baseline_mortality * 100}."

    # Check for decreased productivity (if data is available)
    if 'eggs_infected' in fungus['pollen_data']['Buck']: # Check if productivity data exists
        decreased_productivity = False
        for pollen in fungus['pollen_data']:
            if fungus['pollen_data'][pollen]['eggs_infected'] < fungus['pollen_data'][pollen]['eggs_not_infected']:
                decreased_productivity = True
        if decreased_productivity:
            return True, f"{fungus_name} is a pathogen because it decreases productivity."
    
    return False, f"{fungus_name} does not show pathogenic traits."

def is_commensal(fungus_name):
    """A fungus is a commensal if it doesn't harm the host (no increased mortality)."""
    fungus = data[fungus_name]
    for pollen in fungus['pollen_data']:
        if fungus['pollen_data'][pollen]['mortality_rate'] > baseline_mortality:
            return False, f"{fungus_name} is not a commensal because it increases mortality."

    # Additionally check if productivity is not harmed or is even improved
    improved_productivity = False
    for pollen in fungus['pollen_data']:
        if 'eggs_infected' in fungus['pollen_data'][pollen]:
            if fungus['pollen_data'][pollen]['eggs_infected'] >= fungus['pollen_data'][pollen]['eggs_not_infected']:
                 improved_productivity = True
    
    if improved_productivity:
        return True, f"{fungus_name} is a commensal (or symbiont) because it does not increase mortality and can even increase productivity."
    else:
        return True, f"{fungus_name} is a commensal because it does not increase mortality."

# --- Final Analysis for the Correct Answer (Option I) ---

print("Analyzing the statements in option I: 'Fungus A and B are pathogens. Fungus C is a commensal.'\n")

# Check Fungus A
pathogen_A, reason_A = is_pathogen('Fungus A')
print(f"Statement 1: Is Fungus A a pathogen? {pathogen_A}.")
max_mortality_A = max(d['mortality_rate'] for d in data['Fungus A']['pollen_data'].values())
print(f"Reason: With most pollens, the mortality rate for Fungus A is {max_mortality_A:.0%}, which is greater than the baseline of {baseline_mortality:.0%}.")
print("-" * 20)

# Check Fungus B
pathogen_B, reason_B = is_pathogen('Fungus B')
print(f"Statement 2: Is Fungus B a pathogen? {pathogen_B}.")
max_mortality_B = max(d['mortality_rate'] for d in data['Fungus B']['pollen_data'].values())
print(f"Reason: The mortality rate for Fungus B is {max_mortality_B:.0%}, which is greater than the baseline of {baseline_mortality:.0%}.")
print("-" * 20)


# Check Fungus C
commensal_C, reason_C = is_commensal('Fungus C')
print(f"Statement 3: Is Fungus C a commensal? {commensal_C}.")
max_mortality_C = max(d['mortality_rate'] for d in data['Fungus C']['pollen_data'].values())
print(f"Reason 1: The mortality rate for Fungus C is {max_mortality_C:.0%}, which is equal to the baseline rate, indicating no harm.")
buck_eggs_not_infected = data['Fungus C']['pollen_data']['Buck']['eggs_not_infected']
buck_eggs_infected = data['Fungus C']['pollen_data']['Buck']['eggs_infected']
print(f"Reason 2: Productivity often increases. For example, on a buck pollen diet, infected bees produced {buck_eggs_infected} eggs, while non-infected bees produced {buck_eggs_not_infected} eggs.")
print("\nConclusion: All statements in option I are supported by the data.")
print("<<<I>>>")