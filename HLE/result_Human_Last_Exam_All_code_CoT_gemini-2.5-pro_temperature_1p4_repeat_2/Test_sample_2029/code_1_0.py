import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), or n-choose-k.
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_match_sequences():
    """
    Calculates the total number of different score sequences (a_1,...,a_8, b_1,...,b_6)
    by finding the coefficient of x^48 in the expansion of a product of generating functions.
    The formula is: Sum_{i=0..8, j=0..6} (-1)^(i+j) * C(8,i) * C(6,j) * C(61 - 7i - 9j, 13).
    """
    print("Calculating the total number of sequences based on the formula:")
    print("Sum_{i,j} [(-1)^(i+j) * C(8,i) * C(6,j) * C(61 - 7*i - 9*j, 13)]\n")

    total_sequences = 0
    
    # Iterate i from 0 to 8 for the term (1-x^7)^8
    # Iterate j from 0 to 6 for the term (1-x^9)^6
    for i in range(9):
        for j in range(7):
            
            # The power of x contributed by the first two terms is 7*i + 9*j
            # The remaining power needed from (1-x)^-14 is k.
            k = 48 - (7 * i + 9 * j)
            
            # If k is negative, this combination of i and j is not possible for x^48
            if k < 0:
                continue

            # The coefficient from expanding (1-x)^-14 is C(k + 14 - 1, 14 - 1)
            # This is equivalent to C(61 - 7i - 9j, 13)
            n_comb_c = 61 - 7*i - 9*j
            k_comb_c = 13
            
            # C(n,k) is zero if n < k. This check prunes terms that are zero.
            if n_comb_c < k_comb_c:
                continue

            # Calculate the components of the term in the sum
            coeff_a = combinations(8, i)
            coeff_b = combinations(6, j)
            coeff_c = combinations(n_comb_c, k_comb_c)
            sign = (-1)**(i + j)
            
            term_value = sign * coeff_a * coeff_b * coeff_c

            # Only print non-zero terms
            if term_value == 0:
                continue
            
            # Output the numbers used in the calculation for this term
            sign_str = "+" if sign > 0 else "-"
            print(f"Term (i={i}, j={j}): {sign_str} [C(8,{i})*C(6,{j})*C({n_comb_c},{k_comb_c})] = {sign_str} [{coeff_a}*{coeff_b}*{coeff_c}] = {term_value}")

            total_sequences += term_value

    print("\n-------------------------------------------------------------")
    print(f"The final sum of all terms gives the total number of sequences:")
    print(total_sequences)


solve_match_sequences()
<<<2611776918150>>>