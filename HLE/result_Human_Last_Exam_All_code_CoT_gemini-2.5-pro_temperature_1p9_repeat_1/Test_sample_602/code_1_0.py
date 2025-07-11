import math

def calculate_l(n):
    """
    Calculates the exact value of l(n) for n >= 5.

    Args:
        n: An integer greater than or equal to 5.
    
    Returns:
        The numerical value of l(n).
    """
    if n < 5:
        raise ValueError("n must be greater than or equal to 5")

    # Let's denote the two values in the definition of M_(n) as 'a' and 'b'
    # a corresponds to the diagonal elements, b to the off-diagonal
    a_sq = 1 - (n - 1) / (n**2)
    a = math.sqrt(a_sq)
    b = 1 / n
    
    # We derived the expression for l(n) as:
    # l(n) = 6 - 2 * (2*a^2 + (2*n-1)*a*b + (2*n-3)*b^2)
    # Let's calculate the term in the parenthesis
    
    term_a_sq = 2 * a_sq
    term_ab = (2 * n - 1) * a * b
    term_b_sq = (2 * n - 3) * b**2
    
    total_term = term_a_sq + term_ab + term_b_sq
    
    l_n = 6 - 2 * total_term
    
    # Printing the components of the formula
    print(f"For n = {n}:")
    print(f"  a = sqrt(1 - (n-1)/n^2) = {a}")
    print(f"  b = 1/n = {b}")
    print("\n  Components of the expression T = 2*a^2 + (2*n-1)*a*b + (2*n-3)*b^2:")
    print(f"  2*a^2 = {term_a_sq}")
    print(f"  (2*n-1)*a*b = {term_ab}")
    print(f"  (2*n-3)*b^2 = {term_b_sq}")
    print(f"  T = {total_term}")
    
    print("\n  Final calculation l(n) = 6 - 2*T:")
    print(f"  l({n}) = {l_n}")
    
    # As a cross-check, we use the simplified analytical formula:
    # l(n) = 2 + 2/n^2 - (2*(2*n-1)*sqrt(n^2-n+1))/n^2
    l_n_simplified = 2 + 2/n**2 - (2 * (2*n - 1) * math.sqrt(n**2 - n + 1)) / n**2
    print(f"\nValue from simplified formula: {l_n_simplified}")
    
    return l_n

if __name__ == '__main__':
    # We are asked to calculate the value for n >= 5. Let's use n=5.
    n_val = 5
    final_value = calculate_l(n_val)
