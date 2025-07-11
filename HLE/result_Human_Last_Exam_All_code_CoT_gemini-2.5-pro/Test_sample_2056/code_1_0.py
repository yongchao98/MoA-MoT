import sympy as sp

def solve_lk_n():
    """
    This function defines the derived formula for l_k(n) using sympy and prints it.
    The derivation involves several steps from hyperbolic geometry and linear algebra.
    The final expression is a combination of terms involving n and k.
    """
    
    # Define the symbols 'n' and 'k' for our symbolic computation.
    # n is a positive integer (n>=3), and k is a positive real constant.
    n = sp.Symbol('n', integer=True, positive=True)
    k = sp.Symbol('k', real=True, positive=True)

    # The formula for l_k(n) has been derived as follows:
    # l_k(n) = (1/2)*ln(n+1) + (n-1)*ln(k) - k^2 * (2n-1) / (2n)
    # The numbers in this equation are:
    # coefficient for ln(n+1): 1/2
    # constant added to n in ln: 1
    # coefficient for ln(k): (n-1) -> consists of 1 and -1
    # coefficient for k^2 term: -(2n-1)/(2n) -> consists of 2, -1, 2
    
    # Let's construct the expression term by term.
    
    # Term 1: from the determinant of the covariance matrix Sigma
    # This term is (1/2) * ln(n + 1)
    term1_coeff = sp.Rational(1, 2)
    term1 = term1_coeff * sp.log(n + 1)

    # Term 2: from the volume correction factor in the hyperbolic PDF
    # This term is (n - 1) * ln(k)
    term2 = (n - 1) * sp.log(k)

    # Term 3: from the main exponential part of the Normal distribution PDF
    # This term is -k^2 * (2*n - 1) / (2*n)
    term3_coeff_num = 2 * n - 1
    term3_coeff_den = 2 * n
    term3 = - (k**2 * term3_coeff_num) / term3_coeff_den
    
    # The final expression for l_k(n) is the sum of these three terms.
    l_k_n_expression = term1 + term2 + term3

    # Print the final result in a formatted way.
    print("The exact value of l_k(n) is given by the expression:")
    
    # We can use sympy's pretty print for a nice mathematical layout.
    sp.pprint(l_k_n_expression, use_unicode=True)

if __name__ == '__main__':
    solve_lk_n()