def solve_annihilation_time(N1, M1, N2, M2):
    """
    Calculates the expected time for the second collision of five particles.

    Args:
        N1, M1, N2, M2: Positive integers defining the initial particle gaps.
    """

    # The initial gaps between adjacent particles are:
    # d1 = N1
    # d2 = M1
    # d3 = N2
    # d4 = M2

    # The expectation of the time tau is given by the formula:
    # E[tau] = 1/2 * (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
    
    term1 = N1 * M1
    term2 = N1 * N2
    term3 = N1 * M2
    term4 = M1 * N2
    term5 = M1 * M2
    term6 = N2 * M2
    
    sum_of_products = term1 + term2 + term3 + term4 + term5 + term6
    
    expected_time = 0.5 * sum_of_products
    
    # The problem asks to output the final equation, showing each number.
    print(f"The initial gaps are N1={N1}, M1={M1}, N2={N2}, M2={M2}.")
    print("The expected time E[tau] is calculated by the formula:")
    print("E[tau] = 1/2 * (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2)")
    print(f"E[tau] = 1/2 * ({N1}*{M1} + {N1}*{N2} + {N1}*{M2} + {M1}*{N2} + {M1}*{M2} + {N2}*{M2})")
    print(f"E[tau] = 1/2 * ({term1} + {term2} + {term3} + {term4} + {term5} + {term6})")
    print(f"E[tau] = 1/2 * ({sum_of_products})")
    print(f"E[tau] = {expected_time}")

# Example usage with some placeholder values for N1, M1, N2, M2.
# The user should replace these with their specific values.
if __name__ == '__main__':
    # As an example, let's use N1=1, M1=2, N2=3, M2=4
    N1 = 1
    M1 = 2
    N2 = 3
    M2 = 4
    solve_annihilation_time(N1, M1, N2, M2)
