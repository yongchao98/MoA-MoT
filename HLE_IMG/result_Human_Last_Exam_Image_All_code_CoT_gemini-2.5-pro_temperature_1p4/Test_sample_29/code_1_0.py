import numpy as np

def calculate_hill_numbers():
    """
    Calculates and prints Hill numbers for a given set of species counts.
    """
    # Species counts identified from the image
    counts = np.array([62, 22, 6])
    q_values = [1, 2, 3, 4]

    # Total number of individuals
    N = np.sum(counts)
    # Proportional abundances
    p = counts / N
    
    results = []

    print("--- Hill Number Calculations ---")
    print(f"Species counts (n_i): {counts.tolist()}")
    print(f"Total individuals (N): {N}\n")

    # Calculate for q=1 (Shannon Diversity)
    q = 1
    shannon_entropy = -np.sum(p * np.log(p))
    hill_1 = np.exp(shannon_entropy)
    results.append(hill_1)
    
    # Print the equation for q=1
    sum_terms_str = ' + '.join([f"({ni}/{N})*ln({ni}/{N})" for ni in counts])
    print(f"For q={q}:")
    print(f"D^{q} = exp(-[ {sum_terms_str} ])")
    print(f"D^{q} = {hill_1:.2f}\n")

    # Calculate for q > 1
    for q in q_values[1:]:
        sum_p_q = np.sum(p**q)
        hill_q = (sum_p_q)**(1 / (1 - q))
        results.append(hill_q)

        # Print the equation for q
        sum_terms_str = ' + '.join([f"({ni}/{N})^{q}" for ni in counts])
        print(f"For q={q}:")
        print(f"D^{q} = ( {sum_terms_str} )^(1/(1-{q}))")
        print(f"D^{q} = {hill_q:.2f}\n")

    # Print the final comma-separated result
    formatted_results = [f"{res:.2f}" for res in results]
    print("--- Final Result ---")
    print("Hill numbers for q=1, 2, 3, 4:")
    print(",".join(formatted_results))

# Execute the calculation
calculate_hill_numbers()