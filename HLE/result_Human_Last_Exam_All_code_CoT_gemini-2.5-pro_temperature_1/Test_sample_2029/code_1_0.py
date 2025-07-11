import math

def combinations(n, k):
    """
    Calculates the number of combinations (n choose k).
    Returns 0 if k > n or k < 0, as is standard in combinatorics.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the total number of different sequences (a_1,...,a_8, b_1,...,b_6)
    based on the derived formula from generating functions.
    """
    total_count = 0
    
    print("The final calculation is based on the formula:")
    print("Sum_{i=0 to 8} Sum_{j=0 to 6} C(8,i)*C(6,j)*(-1)^(i+j)*C(61-7i-9j, 13)\n")
    print("Calculating the non-zero terms of the sum:\n")

    # Iterate through the first binomial expansion (1-x^7)^8
    for i in range(9):  # i from 0 to 8
        # Iterate through the second binomial expansion (1-x^9)^6
        for j in range(7):  # j from 0 to 6
            
            # The power of x from these two expansions
            power_from_ij = 7 * i + 9 * j
            
            # We need the coefficient from (1-x)^-14 for x^(48 - power_from_ij)
            # which is C( (48-power_from_ij) + 13, 13 )
            n_comb = 61 - power_from_ij
            k_comb = 13
            
            # Calculate the components of the term
            term_c8i = combinations(8, i)
            term_c6j = combinations(6, j)
            sign = (-1)**(i + j)
            term_cnk = combinations(n_comb, k_comb)
            
            # The full term in the summation
            term_value = term_c8i * term_c6j * sign * term_cnk
            
            if term_value != 0:
                # Print the calculation for the current term
                sign_str = "+" if sign > 0 else "-"
                print(f"Term (i={i}, j={j}): {sign_str} C(8,{i})*C(6,{j})*C({n_comb},{k_comb}) = {term_value}")
                total_count += term_value

    print(f"\nThe total number of different sequences is the sum of these terms.")
    print(f"\nFinal Answer: {total_count}")

solve()