import numpy as np

def solve_system():
    """
    This function calculates the expected number of rounds R and provides a detailed breakdown of the calculation.
    """
    # The system's evolution is a linear transformation over F_2, S_new = T * S_old.
    # The matrix T has been derived in the explanation. T = I + N.
    # The order of T is 8, so cycle lengths can be 1, 2, 4, 8.
    
    # Let n(d) be the number of states S such that T^d * S = S.
    # This is equal to 2^dim(Ker(T^d - I)).
    # For d being a power of 2, T^d - I = N^d (mod 2).
    # Based on solving the systems of linear equations N^d * x = 0:
    dim_ker_N1 = 1
    dim_ker_N2 = 2
    dim_ker_N4 = 4
    dim_ker_N8 = 8  # N^8 is the zero matrix

    n_1 = 2**dim_ker_N1
    n_2 = 2**dim_ker_N2
    n_4 = 2**dim_ker_N4
    n_8 = 2**dim_ker_N8

    # Let c(d) be the number of states with a cycle length of exactly d.
    # We can find these by inclusion-exclusion.
    c_1 = n_1
    c_2 = n_2 - n_1
    c_4 = n_4 - n_2
    c_8 = n_8 - n_4
    
    total_states = 2**8
    
    print("The expected number of rounds, E[R], is the average cycle length over all 256 possible initial states.")
    print("The possible cycle lengths are 1, 2, 4, and 8.")
    
    print("\nNumber of states for each exact cycle length:")
    print(f"c(1): {c_1}")
    print(f"c(2): {c_2}")
    print(f"c(4): {c_4}")
    print(f"c(8): {c_8}")

    # E[R] = (Sum of all cycle lengths) / (Total number of states)
    term1 = 1 * c_1
    term2 = 2 * c_2
    term3 = 4 * c_4
    term4 = 8 * c_8
    sum_of_lengths = term1 + term2 + term3 + term4
    
    expected_R = sum_of_lengths / total_states

    print(f"\nThe total sum of all cycle lengths is calculated as follows:")
    print(f"Sum = (1 * {c_1}) + (2 * {c_2}) + (4 * {c_4}) + (8 * {c_8})")
    print(f"Sum = {term1} + {term2} + {term3} + {term4} = {sum_of_lengths}")
    
    print(f"\nThe expected value E[R] is this sum divided by the total number of states ({total_states}):")
    print(f"E[R] = {sum_of_lengths} / {total_states} = {expected_R}")

    final_answer = round(expected_R, 2)
    print(f"\nFinally, rounded to two decimal places, the expected value is: {final_answer}")

solve_system()