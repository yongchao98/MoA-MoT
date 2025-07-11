import math

def calculate_sequences():
    """
    Calculates the number of possible sequences based on user-provided N, K, and M.

    A sequence is defined as K positive integers a_1, a_2, ..., a_K such that:
    1. 1 <= a_1 < a_2 < ... < a_K <= N
    2. a_{i+1} - a_i <= M for i = 1, ..., K-1
    3. The condition M*(K-1) < N holds.
    """
    try:
        N_str = input("Enter the maximum value for an element (N): ")
        K_str = input("Enter the length of the sequence (K): ")
        M_str = input("Enter the maximum increase between consecutive elements (M): ")

        N = int(N_str)
        K = int(K_str)
        M = int(M_str)

        if not (N > 0 and K > 0 and M > 0):
            print("Error: N, K, and M must be positive integers.")
            return

        # As per the problem description, this condition is given to hold.
        # if M * (K - 1) >= N:
        #     print(f"Warning: The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        #     print("The formula may not be directly applicable or may yield a count of 0.")

        total_sequences = 0
        equation_symbolic = []
        equation_numeric = []
        term_values = []

        # The calculation is based on the inclusion-exclusion principle:
        # Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
        for j in range(K):
            n_for_comb2 = N - j * M
            
            # math.comb(n, k) is 0 if k > n, which correctly handles terms where n_for_comb2 < K.
            c1_val = math.comb(K - 1, j)
            c2_val = math.comb(n_for_comb2, K)
            
            term_val = c1_val * c2_val
            
            # Symbolic representation of the term
            c1_str = f"C({K-1}, {j})"
            c2_str = f"C({n_for_comb2}, {K})"
            symbolic_part = f"{c1_str} * {c2_str}"
            
            # Numeric representation of the term
            numeric_part = f"{c1_val} * {c2_val}"

            if j % 2 == 1: # Odd j corresponds to subtraction
                total_sequences -= term_val
                equation_symbolic.append(f"- {symbolic_part}")
                equation_numeric.append(f"- {numeric_part}")
            else: # Even j corresponds to addition
                total_sequences += term_val
                if j > 0:
                    equation_symbolic.append(f"+ {symbolic_part}")
                    equation_numeric.append(f"+ {numeric_part}")
                else: # First term has no preceding sign
                    equation_symbolic.append(symbolic_part)
                    equation_numeric.append(numeric_part)

        print("\nThe number of possible sequences is calculated using the formula derived from the Principle of Inclusion-Exclusion:")
        print("Number of sequences = " + " ".join(equation_symbolic))
        print("                    = " + " ".join(equation_numeric))
        print(f"                    = {total_sequences}")

    except ValueError:
        print("Invalid input. Please enter valid integers for N, K, and M.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_sequences()