import math

def calculate_destabilizers(n):
    """
    Calculates the number of destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n} and prints the detailed calculation.
    
    The formula is: 4^n * 2^(n*(n+1)/2)
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Using integer arithmetic to avoid float precision issues with large numbers
    
    # Part 1 of the formula: 4^n
    term1_base = 4
    # Use pow for large integer exponentiation
    term1_val = pow(term1_base, n)

    # Part 2 of the formula: 2^(n*(n+1)/2)
    term2_base = 2
    exp_n_plus_1 = n + 1
    exp_numerator = n * exp_n_plus_1
    # Use integer division
    exp_val = exp_numerator // 2
    # Use pow for large integer exponentiation
    term2_val = pow(term2_base, exp_val)

    # Final result
    total = term1_val * term2_val

    # Print the derivation step-by-step
    print(f"To find the number of destabilizer sets for n = {n}:")
    print(f"The formula is: {term1_base}^n * {term2_base}^(n * (n + 1) / 2)")
    print("Substituting n with the given value:")
    print(f"= {term1_base}^{n} * {term2_base}^({n} * ({n} + 1) / 2)")
    print(f"= {term1_base}^{n} * {term2_base}^({n} * {exp_n_plus_1} / 2)")
    print(f"= {term1_val} * {term2_base}^({exp_numerator} / 2)")
    print(f"= {term1_val} * {term2_base}^{exp_val}")
    print(f"= {term1_val} * {term2_val}")
    print("---------------------------------")
    print(f"The total number of destabilizer sets is: {total}")

# --- User-configurable value for n ---
# You can change this value to calculate for a different number of qubits.
n_qubits = 4
# ---

calculate_destabilizers(n_qubits)
