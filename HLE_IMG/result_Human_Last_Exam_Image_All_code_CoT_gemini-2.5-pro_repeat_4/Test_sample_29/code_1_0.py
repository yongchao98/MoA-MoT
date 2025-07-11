import numpy as np

def calculate_and_print_hill(q, abundances, name):
    """
    Calculates the Hill number for a given q and list of abundances,
    and prints the steps of the calculation.
    """
    print(f"--- Calculation for Taxonomic Level: {name} (q={q}) ---")
    
    # Convert to numpy array for vectorized operations
    abundances = np.array(abundances)
    total_individuals = np.sum(abundances)
    
    # Calculate proportional abundances, filtering out any groups with zero individuals
    proportions = abundances[abundances > 0] / total_individuals
    
    print(f"Abundance counts (n_i): {abundances}")
    print(f"Total individuals (N): {total_individuals}")
    print("Proportional abundances (p_i): [" + ", ".join([f"{p:.4f}" for p in proportions]) + "]")

    # Handle the special case for q = 1
    if q == 1:
        # D_1 = exp(-Σ(p_i * ln(p_i)))
        shannon_entropy = -np.sum(proportions * np.log(proportions))
        hill_number = np.exp(shannon_entropy)
        
        print("Equation: D_1 = exp(-Σ(p_i * ln(p_i)))")
        print(f"Value of Σ(p_i * ln(p_i)): {-shannon_entropy:.4f}")
        print(f"Hill Number (D_1): exp({shannon_entropy:.4f}) = {hill_number:.4f}")

    # Handle the general case for q != 1
    else:
        # D_q = (Σ(p_i^q))^(1/(1-q))
        sum_of_proportions_raised_to_q = np.sum(proportions**q)
        exponent = 1 / (1 - q)
        hill_number = sum_of_proportions_raised_to_q**exponent
        
        print(f"Equation: D_{q} = (Σ(p_i^{q}))^(1 / (1 - {q}))")
        sum_eq_str = " + ".join([f"{p:.4f}^{q}" for p in proportions])
        print(f"Value of Σ(p_i^{q}) = {sum_eq_str} = {sum_of_proportions_raised_to_q:.4f}")
        print(f"Hill Number (D_{q}): ({sum_of_proportions_raised_to_q:.4f})^({exponent:.2f}) = {hill_number:.4f}")
    
    print("-" * (len(name) + 33))
    print() # Add a blank line for readability
    
    return hill_number

# Step 1: Define abundances based on visual identification and counting
abundances_order = [97]
abundances_family = [97]
abundances_genus = [97]
abundances_species = [64, 22, 6, 5]

# Step 2: Define parameters for each calculation
calculation_params = [
    {'name': 'Order',   'q': 1, 'abundances': abundances_order},
    {'name': 'Family',  'q': 2, 'abundances': abundances_family},
    {'name': 'Genus',   'q': 3, 'abundances': abundances_genus},
    {'name': 'Species', 'q': 4, 'abundances': abundances_species}
]

# Step 3: Perform calculations and store the results
results = []
for params in calculation_params:
    result = calculate_and_print_hill(params['q'], params['abundances'], params['name'])
    results.append(result)

# Step 4: Format and print the final answer
rounded_results = [f"{res:.2f}" for res in results]
final_answer_string = ", ".join(rounded_results)

print("===============================================")
print("Final Hill Numbers (Order, Family, Genus, Species):")
print(final_answer_string)
print("===============================================")
<<<1.00, 1.00, 1.00, 1.73>>>