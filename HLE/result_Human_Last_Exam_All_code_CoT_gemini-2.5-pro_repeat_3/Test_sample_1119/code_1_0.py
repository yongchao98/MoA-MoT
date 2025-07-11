import math

def count_sequences(N, K, M):
    """
    Calculates the number of sequences of K positive integers up to N,
    where each number is greater than the one before it, and the increase
    between consecutive numbers does not exceed M.

    The calculation is based on the formula derived from inclusion-exclusion:
    Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - M*j, K)]
    """
    # Verify the condition given in the problem statement
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print(f"The result is calculated, but the problem's premise is violated.\n")

    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}")
    print(f"The formula used is: Sum_{{j=0 to {K-1}}} [(-1)^j * C({K-1}, j) * C({N} - {M}*j, {K})]")
    print("-" * 30)

    total_sequences = 0
    final_equation = []

    for j in range(K):
        # Calculate the first binomial coefficient: C(K-1, j)
        comb1 = math.comb(K - 1, j)

        # Calculate the second binomial coefficient: C(N - M*j, K)
        n_val = N - M * j
        k_val = K

        # math.comb(n, k) returns 0 if k > n, which is the desired behavior.
        # It raises a ValueError for negative n, so we handle that case.
        if n_val < 0:
            comb2 = 0
        else:
            comb2 = math.comb(n_val, k_val)

        # Calculate the term for the current j
        term = ((-1)**j) * comb1 * comb2

        # Build the string for the detailed calculation step
        sign_char = "-" if j % 2 != 0 else "+"
        term_str = (
            f"j={j}: (-1)^{j} * C({K-1}, {j}) * C({N}-{M}*{j}, {K}) "
            f"= {int((-1)**j)} * {comb1} * {comb2} = {term}"
        )
        print(term_str)

        # Append the term to the final equation string and update the total
        if j == 0:
            final_equation.append(f"{term}")
        else:
            final_equation.append(f"{sign_char} {abs(term)}")
        
        total_sequences += term

    print("-" * 30)
    print(f"Final Calculation: {' '.join(final_equation)} = {total_sequences}")
    print(f"The total number of possible sequences is: {total_sequences}")
    return total_sequences

if __name__ == '__main__':
    # User-provided parameters from the problem description
    # As no specific values were given, we use an example where M*(K-1) < N holds.
    N_val = 10
    K_val = 4
    M_val = 3
    final_answer = count_sequences(N_val, K_val, M_val)
    print(f"\n<<< {final_answer} >>>")
