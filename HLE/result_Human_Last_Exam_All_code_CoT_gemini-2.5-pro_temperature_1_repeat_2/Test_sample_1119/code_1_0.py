import math

def calculate_sequences():
    """
    Calculates the number of possible sequences based on user-provided N, K, and M.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. Each number is at most N.
    2. Each number is greater than the one before it (a_i > a_{i-1}).
    3. The increase between consecutive numbers does not exceed M (a_{i+1} - a_i <= M).
    4. A precondition M*(K-1) < N holds.
    """
    try:
        N = int(input("Enter the maximum value for any number in the sequence (N): "))
        K = int(input("Enter the length of the sequence (K): "))
        M = int(input("Enter the maximum increase between consecutive numbers (M): "))

        if not (N > 0 and K > 0 and M > 0):
            print("\nError: N, K, and M must be positive integers.")
            return

        # The problem states the condition M*(K-1) < N holds.
        # We can check it to ensure the inputs are as expected.
        if M * (K - 1) >= N:
            print(f"\nWarning: The condition M*(K-1) < N is not met for your inputs.")
            print(f"M*(K-1) = {M*(K-1)}, which is not less than N = {N}.")
            print("The calculation will proceed, but the result might not be meaningful under the problem's specified conditions.")

        # Helper for combinations C(n, r)
        def nCr(n, r):
            if r < 0 or r > n:
                return 0
            return math.comb(n, r)

        total_sequences = 0
        symbolic_parts = []
        numeric_parts = []
        
        # Formula: Sum_{p=0 to K-1} [(-1)^p * C(K-1, p) * C(N - p*M, K)]
        for p in range(K):
            sign = (-1)**p
            comb1_n, comb1_r = K - 1, p
            comb2_n, comb2_r = N - p * M, K

            comb1_val = nCr(comb1_n, comb1_r)
            comb2_val = nCr(comb2_n, comb2_r)

            # Skip terms that are zero to keep the equation clean
            if comb1_val == 0 or comb2_val == 0:
                continue

            term_value = sign * comb1_val * comb2_val
            total_sequences += term_value

            # Build the string for the equation
            op = " + " if sign > 0 else " - "
            
            symbolic_part = f"C({comb1_n}, {comb1_r}) * C({comb2_n}, {comb2_r})"
            numeric_part = f"{comb1_val} * {comb2_val}"

            if not symbolic_parts: # First term
                 symbolic_parts.append(symbolic_part)
                 numeric_parts.append(numeric_part)
            else:
                 symbolic_parts.append(op + symbolic_part)
                 numeric_parts.append(op + numeric_part)


        print("\n----------------------------------------------------")
        print("The number of sequences is calculated using the formula:")
        print("Count = Sum_{p=0 to K-1} [(-1)^p * C(K-1, p) * C(N - p*M, K)]")
        print("\nFor your inputs N={}, K={}, M={}:".format(N, K, M))
        
        if not symbolic_parts:
            # This happens if the first term itself is zero, e.g., N < K
            final_equation = "0"
            final_numeric = "0"
        else:
            final_equation = "".join(symbolic_parts)
            final_numeric = "".join(numeric_parts)
            # handle leading " + " for the first term if p=0 was skipped
            if final_equation.startswith(" + "):
                final_equation = final_equation[3:]
            if final_numeric.startswith(" + "):
                final_numeric = final_numeric[3:]

        print(f"\nEquation:")
        print(final_equation)
        print(f"\n= {final_numeric}")
        print(f"\n= {total_sequences}")
        print("----------------------------------------------------")
        
    except ValueError:
        print("\nError: Please enter valid integers for N, K, and M.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Run the calculation
calculate_sequences()