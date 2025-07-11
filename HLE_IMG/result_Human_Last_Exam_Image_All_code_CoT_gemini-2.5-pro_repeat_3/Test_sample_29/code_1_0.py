import numpy as np

def calculate_and_print_hill_number(abundances, q):
    """Calculates and prints the Hill number of order q for given abundances."""
    print(f"--- Calculating Hill Number for q = {q} ---")
    N = np.sum(abundances)
    p = abundances / N

    print(f"Species Abundances (n_i): {abundances}")
    print(f"Total Individuals (N): {N}")
    print(f"Proportions (p_i): {[f'{pi:.4f}' for pi in p]}")

    if q == 1:
        # Special case for q=1: D_1 = exp(-sum(p_i * ln(p_i)))
        # Each number in the equation:
        p_log_p_terms = p * np.log(p)
        shannon_entropy_sum = np.sum(p_log_p_terms)
        hill_number = np.exp(-shannon_entropy_sum)
        
        print("\nEquation: D_1 = exp(- sum(p_i * ln(p_i)))")
        print("Individual p_i*ln(p_i) terms: " + ", ".join([f"{term:.4f}" for term in p_log_p_terms]))
        print(f"The sum is: {shannon_entropy_sum:.4f}")
        print(f"Final Calculation: D_1 = exp(-({shannon_entropy_sum:.4f}))")
        
    else:
        # General formula for q != 1: D_q = (sum(p_i^q))^(1/(1-q))
        # Each number in the equation:
        p_raised_to_q = p**q
        sum_p_q = np.sum(p_raised_to_q)
        exponent = 1 / (1 - q)
        hill_number = sum_p_q**exponent
        
        print(f"\nEquation: D_{q} = (sum(p_i^{q}))^(1/(1-{q}))")
        print(f"Individual p_i^{q} terms: " + ", ".join([f"{term:.6f}" for term in p_raised_to_q]))
        print(f"The sum is: {sum_p_q:.6f}")
        print(f"Final Calculation: D_{q} = ({sum_p_q:.6f})^({exponent:.2f})")

    print(f"Result for D_{q}: {hill_number:.4f}")
    print("-" * (35 + len(str(q))))
    return hill_number

# Species counts from the image
abundances = np.array([64, 24, 4, 4])

# Orders of q to be calculated
q_values = [1, 2, 3, 4]

# Calculate and store the results
results = []
for q in q_values:
    result = calculate_and_print_hill_number(abundances, q)
    results.append(result)

# Format and print the final consolidated answer
formatted_results = [f"{r:.2f}" for r in results]
final_answer_string = ", ".join(formatted_results)

print("\nFinal consolidated answer:")
print(final_answer_string)
<<<2.41, 1.96, 1.79, 1.71>>>