import numpy as np

def get_continued_fraction_coeffs(p, q):
    """Computes the continued fraction coefficients for a rational number p/q."""
    coeffs = []
    if q == 0:
        return []
    while q != 0:
        coeffs.append(p // q)
        p, q = q, p % q
    return coeffs

def compute_generalized_markov_number(p, q):
    """Computes the generalized Markov number m_{p/q}."""
    if q == 0 or p < 0 or q < 0 or p >= q:
        raise ValueError("Input rational must be in the interval [0, 1)")
    
    # Step 1: Get continued fraction coefficients for p/q
    cf_coeffs = get_continued_fraction_coeffs(p, q)
    # We need the coefficients a_1, a_2, ... from [0; a_1, a_2, ...]
    if not cf_coeffs or cf_coeffs[0] != 0:
         raise ValueError("The rational number should be in [0, 1) for this definition.")
    a = cf_coeffs[1:]
    
    # Step 2: Define matrices L and R
    L = np.array([[2, 1], [1, 1]])
    R = np.array([[1, 1], [1, 2]])
    
    # Step 3: Construct the matrix M
    M = np.identity(2)
    matrices = [L, R]
    
    for i, power in enumerate(a):
        # Alternate between L and R
        matrix_to_use = matrices[i % 2]
        M = M @ np.linalg.matrix_power(matrix_to_use, power)
        
    # Step 4: Calculate the Markov number
    trace = np.trace(M)
    # The trace should be divisible by 3
    if trace % 3 != 0:
        # This case should ideally not happen for valid p/q
        raise ValueError("Trace is not divisible by 3")
        
    markov_number = int(trace / 3)
    return markov_number

# Main execution
# Compute the generalized Markov number for p/q = 4/7
p, q = 4, 7
m_4_7 = compute_generalized_markov_number(p, q)

# Compute the continued fraction of the resulting number m_4_7
final_cf_coeffs = get_continued_fraction_coeffs(m_4_7, 1)

# Print the result
# The problem asks to output each number in the final equation.
# The continued fraction of 34 is [34]. We will print the coefficients.
print(f"The generalized Markov number m_{{{p}/{q}}} is: {m_4_7}")
print(f"The continued fraction of {m_4_7} is {final_cf_coeffs}.")
print("The numbers in the final equation (the coefficients of the continued fraction) are:")
for coeff in final_cf_coeffs:
    print(coeff)
