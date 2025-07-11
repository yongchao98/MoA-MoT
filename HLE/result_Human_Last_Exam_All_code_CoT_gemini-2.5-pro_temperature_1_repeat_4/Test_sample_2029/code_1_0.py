import math

def combinations(n, k):
    """Computes the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_match_sequences():
    """
    Calculates the number of possible score sequences based on the generating function method.
    """
    total_sequences = 0
    equation_terms = []
    
    # The number of sequences is the coefficient of x^48 in the expansion of
    # (1+...+x^6)^8 * (1+...+x^8)^6.
    # The formula for this coefficient is:
    # Sum_{i=0 to 8} Sum_{j=0 to 6} (-1)^(i+j) * C(8,i) * C(6,j) * C(61 - 7i - 9j, 13)
    
    for i in range(8 + 1):
        for j in range(6 + 1):
            # We only need to compute for terms where the combination is valid
            if 7 * i + 9 * j <= 48:
                c_8_i = combinations(8, i)
                c_6_j = combinations(6, j)
                c_N_K = combinations(61 - 7 * i - 9 * j, 13)
                
                term_value = ((-1)**(i + j)) * c_8_i * c_6_j * c_N_K
                
                if term_value != 0:
                    total_sequences += term_value
                    
                    sign = "-" if term_value < 0 else "+"
                    # Format the term string like: +/- C(n1,k1)*C(n2,k2)*C(n3,k3)
                    term_str = f" {sign} C(8,{i})*C(6,{j})*C({61 - 7*i - 9*j},13)"
                    equation_terms.append(term_str)

    print("The total number of sequences is given by the sum:")
    # Print the equation with all its terms
    # The first term's sign is handled by removing the leading space and '+'
    full_equation = "".join(equation_terms).lstrip(" +")
    print(full_equation)
    
    print("\nWhich evaluates to:")
    print(total_sequences)


solve_match_sequences()