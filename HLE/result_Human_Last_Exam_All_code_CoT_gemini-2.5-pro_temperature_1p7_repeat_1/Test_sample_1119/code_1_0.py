import math

def count_sequences(N, K, M):
    """
    Calculates the number of sequences of K positive integers up to N, where each number is greater than the one before it, and the increase between consecutive numbers does not exceed M.

    The function prints the derivation based on the formula:
    Count = Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    
    Args:
        N (int): The maximum value for any number in the sequence.
        K (int): The length of the sequence.
        M (int): The maximum allowed increase between consecutive numbers.
    """
    if not (K > 0 and N > 0 and M > 0):
        print("N, K, and M must be positive integers.")
        return

    total_count = 0
    equation_parts = []
    calculation_parts = []

    for j in range(K):  # Summation over j from 0 to K-1
        sign = (-1)**j
        
        # Calculate the first combination: C(K-1, j)
        comb1 = math.comb(K - 1, j)
        
        # Calculate the second combination: C(N - j*M, K)
        # The upper value for this combination is N - j*M
        n_for_comb2 = N - j * M
        comb2 = math.comb(n_for_comb2, K)

        # Update the total count
        total_count += sign * (comb1 * comb2)
        
        # For printing: determine the operator sign
        op_symbol = ""
        if j > 0:
            op_symbol = "+ " if sign > 0 else "- "
        
        # Append the parts for the final equation string
        equation_parts.append(f"{op_symbol}C({K-1}, {j}) * C({N - j * M}, {K})")
        calculation_parts.append(f"{op_symbol}{comb1} * {comb2}")
    
    # Print the full explanation and result
    print("The final calculation is based on the Principle of Inclusion-Exclusion.")
    print("Formula: Count = Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]")
    print("-" * 50)
    
    print("Final Equation using the provided values:")
    print("".join(equation_parts))
    
    print("\nCalculation with the numbers plugged in:")
    print("".join(calculation_parts))
    
    print("\nTotal number of possible sequences:")
    print(total_count)


# --- User Input ---
# Please set the values for N, K, and M for your specific problem.
N = 10
K = 4
M = 3

# --- Execution ---
# The code will now calculate the answer based on the values above.
count_sequences(N, K, M)
<<<220>>>