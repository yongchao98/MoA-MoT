import math

def nCr_safe(n, r):
    """
    Calculates the binomial coefficient C(n, r), also known as "n choose r".
    It's defined to be 0 if n < r or r < 0, which is standard in such formulas.
    """
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def solve_sequence_count(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints.

    The number of sequences is given by the formula derived from the principle
    of inclusion-exclusion:
    Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - M*j, K)
    """
    # The problem statement gives this condition to ensure a non-trivial problem.
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")

    total_count = 0
    equation_str_parts = []
    calculation_str_parts = []
    
    # The summation is from j=0 to K-1
    for j in range(K):
        sign = (-1)**j
        
        # Binomial coefficients for the formula
        comb1 = nCr_safe(K - 1, j)
        comb2_n = N - M * j
        comb2 = nCr_safe(comb2_n, K)
        
        # Calculate the value of the current term
        term_value = comb1 * comb2
        
        # Add the signed term to the total
        total_count += sign * term_value
        
        # Build the strings for the equation and calculation steps
        sign_str = " - " if sign < 0 else " + "
        if j == 0:
            sign_str = ""
            
        equation_str_parts.append(f"{sign_str}C({K-1}, {j}) * C({comb2_n}, {K})")
        calculation_str_parts.append(f"{sign_str}{comb1} * {comb2}")
        
    # Assemble the final output strings
    final_equation_str = "Equation:" + "".join(equation_str_parts)
    final_calculation_str = "Calculation:" + "".join(calculation_str_parts)
    
    print(final_equation_str)
    print(f"{final_calculation_str} = {total_count}")
    print(f"\nTotal number of possible sequences: {total_count}")


# --- Set your parameters here ---
# N: The maximum value for any number in the sequence.
# K: The length of the sequence.
# M: The maximum increase between consecutive numbers.

# Example values:
N = 10
K = 4
M = 3

# Run the calculation
solve_sequence_count(N, K, M)
<<<70>>>