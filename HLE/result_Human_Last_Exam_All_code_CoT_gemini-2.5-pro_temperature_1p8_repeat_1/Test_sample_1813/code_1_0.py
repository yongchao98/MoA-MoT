import numpy as np

def get_continued_fraction_coeffs(n, d):
    """
    Computes the regular continued fraction coefficients for n/d.
    Assumes n, d are integers and d > 0.
    """
    coeffs = []
    if d == 0:
        return []
    while d != 0:
        q = n // d
        r = n % d
        coeffs.append(q)
        n, d = d, r
    return coeffs

def compute_generalized_markov_number(k, n):
    """
    Computes the generalized Markov number m_{k/n}.
    This is defined as the trace of a matrix product derived from the
    continued fraction of n/k.
    n/k = [a_1, a_2, ..., a_r]
    m_{k/n} = tr(X_{a_r} * ... * X_{a_1})
    where X_a = [[a, 1], [-1, 0]]
    """
    if k == 0 or n == 0:
        raise ValueError("k and n must be non-zero")
    
    # Get continued fraction coefficients for n/k
    # For k/n = 4/7, we compute for n/k = 7/4
    coeffs = get_continued_fraction_coeffs(n, k)

    # Matrix generator
    def get_X(a):
        return np.array([[a, 1], [-1, 0]], dtype=int)

    # Calculate the matrix product U = X_{a_r} * ... * X_{a_1}
    if not coeffs:
        return None 

    # Initialize U with the last matrix in the product sequence
    # This corresponds to the last coefficient a_r
    U = get_X(coeffs[-1])

    # Multiply from the right by the other matrices in reverse order of coefficients
    for i in range(len(coeffs) - 2, -1, -1):
        U = U @ get_X(coeffs[i])
        
    markov_number = int(np.trace(U))
    return markov_number

def get_cf_for_integer(num):
    """Computes the trivial continued fraction for an integer."""
    return [num]

# Given values
k = 4
n = 7

# Compute the generalized Markov number m_{4/7}
markov_number = compute_generalized_markov_number(k, n)

# Compute the continued fraction of the resulting number
final_cf_coeffs = get_cf_for_integer(markov_number)

# Print the results as requested
print(f"The generalized Markov number m_{k}/{n} is: {markov_number}")
print(f"The continued fraction associated with {markov_number} is: {final_cf_coeffs}")
print(f"The final equation is {markov_number} = {final_cf_coeffs[0]}")
for coeff in final_cf_coeffs:
    print(f"The coefficient in the expansion is: {coeff}")
