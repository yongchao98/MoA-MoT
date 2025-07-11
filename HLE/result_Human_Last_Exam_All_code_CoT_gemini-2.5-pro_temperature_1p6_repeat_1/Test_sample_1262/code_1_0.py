import itertools
from collections import defaultdict

def get_derangement_polynomial_coeffs(n):
    """
    Computes the coefficients of the n-th derangement polynomial d_n(t).

    A derangement is a permutation sigma of {0, 1, ..., n-1} such that sigma[i] != i.
    An excedance is an index i such that sigma[i] > i.
    
    The function returns a dictionary where keys are powers of t (number of excedances)
    and values are the coefficients (number of derangements with that many excedances).
    """
    if n < 2:
        return defaultdict(int)

    coeffs = defaultdict(int)
    # Using 0-based indexing for permutations of {0, 1, ..., n-1}
    elements = range(n)
    for p in itertools.permutations(elements):
        is_derangement = True
        for i in range(n):
            if p[i] == i:
                is_derangement = False
                break
        
        if is_derangement:
            excedances = 0
            for i in range(n):
                if p[i] > i:
                    excedances += 1
            coeffs[excedances] += 1
    return coeffs

def format_polynomial(coeffs):
    """Formats a coefficient dictionary into a string like 'c_k*t^k + ...'"""
    if not coeffs:
        return "0"
    
    # Sort terms by degree in descending order
    sorted_terms = sorted(coeffs.items(), key=lambda item: item[0], reverse=True)
    
    poly_str_parts = []
    for degree, coeff in sorted_terms:
        # Don't print coefficient if it's 1 and degree > 0
        if coeff == 1 and degree > 0:
            coeff_str = ""
        else:
            coeff_str = str(coeff) + ("*" if degree > 0 else "")

        # Don't print t^0 or t^1
        if degree > 1:
            var_str = f"t^{degree}"
        elif degree == 1:
            var_str = "t"
        else:
            var_str = ""
        
        poly_str_parts.append(coeff_str + var_str)
        
    return " + ".join(poly_str_parts)

def analyze_derangement_polynomial(n):
    """Computes and prints the analysis for d_n(t)."""
    print(f"--- Analysis for n={n} ---")
    coeffs = get_derangement_polynomial_coeffs(n)
    
    if not coeffs:
        print("There are no derangements for n < 2.")
        return

    # Polynomial and its properties
    poly_str = format_polynomial(coeffs)
    degree = max(coeffs.keys())
    leading_coeff = coeffs[degree]
    d_n_at_1 = sum(coeffs.values())
    
    print(f"The derangement polynomial d_{n}(t) is: {poly_str}")
    
    # Analysis for part (a)
    degree_lhs = n - 2
    degree_rhs = (n - 1) + degree
    print(f"(a) The degree of H(U_{{{n-1},{n}}})(t) is {n-1}-1 = {degree_lhs}.")
    print(f"    The degree of t^{{{n-1}}}*d_{n}(t) is {n-1} + {degree} = {degree_rhs}.")
    print(f"    Since {degree_lhs} != {degree_rhs}, the identity is false.")
    
    # Analysis for part (b)
    print(f"(b) The degree of d_{n}(t) is {degree}.")
    print(f"    The leading coefficient is {leading_coeff}.")

    # Analysis for part (c) if n=3
    if n == 3:
        print(f"(c) The value of d_{3}(1) is the sum of coefficients: {d_n_at_1}.")
    
    print("-" * 25)

if __name__ == "__main__":
    # Perform analysis for n=3 to answer all parts
    analyze_derangement_polynomial(3)
    
    # Perform analysis for n=4 to further confirm part (b)
    analyze_derangement_polynomial(4)

    # Final concise answer based on the analysis above
    print("\nFinal Answer:")
    print("(a) No; (b) Yes; (c) 2")

    # Final requested format
    final_answer = "(a) No; (b) Yes; (c) 2"
    print(f"\n<<< {final_answer} >>>")