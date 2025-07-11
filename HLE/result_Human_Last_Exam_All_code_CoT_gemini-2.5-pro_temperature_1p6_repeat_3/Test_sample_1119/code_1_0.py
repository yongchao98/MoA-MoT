import math

def solve_sequence_counting():
    """
    Calculates the number of possible sequences based on the given constraints
    using the principle of inclusion-exclusion.
    
    The problem is to find the number of sequences a_1, ..., a_K such that:
    1 <= a_1 < a_2 < ... < a_K <= N
    a_{i+1} - a_i <= M
    
    The derived formula is: Sum_{s=0 to K-1} [(-1)^s * C(K-1, s) * C(N - s*M, K)]
    """
    # Problem parameters
    N = 20
    K = 5
    M = 4

    # Check the given condition M(K-1) < N
    if not M * (K - 1) < N:
        print(f"The condition M(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print(f"Calculation: M*(K-1) = {M* (K-1)}, which is not less than N = {N}.")
        return

    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}\n")

    total_count = 0
    equation_str = "Number of sequences = "
    values_str = ""

    for s in range(K):
        # Binomial coefficient C(K-1, s)
        comb1 = math.comb(K - 1, s)
        
        # Upper part of the second binomial coefficient C(N - s*M, K)
        n_val = N - s * M
        comb2 = math.comb(n_val, K)

        # The value of the current term in the summation
        term_value = comb1 * comb2
        
        # Append to the equation string representation
        term_str = f"C({K-1}, {s})*C({n_val}, {K})"
        
        if s == 0:
            equation_str += term_str
            values_str += str(term_value)
            total_count += term_value
        else:
            if s % 2 == 1: # Odd s, corresponds to subtraction
                equation_str += f" - {term_str}"
                values_str += f" - {term_value}"
                total_count -= term_value
            else: # Even s, corresponds to addition
                equation_str += f" + {term_str}"
                values_str += f" + {term_value}"
                total_count += term_value

    print(equation_str)
    print(f"= {values_str}")
    print(f"= {total_count}")

# Execute the function to solve the problem and print the result.
solve_sequence_counting()