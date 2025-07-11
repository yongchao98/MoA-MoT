import sympy as sp

def solve_lk_n():
    """
    This function calculates and prints the exact value of l_k(n) in terms of n and k.
    """
    # Define n and k as symbolic variables.
    # n is an integer >= 3, k is a constant. We can treat them as generic symbols.
    n, k = sp.symbols('n k')

    # The step-by-step derivation leads to the following expression for l_k(n):
    # l_k(n) = (1/2)*ln(n+1) - k^2 * (2 - 1/n) + (n-1)*ln(k)
    # We construct this expression using sympy.

    # Term 1: from the determinant of Sigma
    term_det = sp.log(n + 1) / 2

    # Term 2: from the quadratic form in the Gaussian exponent
    term_quad = -k**2 * (2 - 1/n)

    # Term 3: from the Jacobian of the exponential map
    term_jac = (n - 1) * sp.log(k)

    # The final expression for l_k(n)
    l_k_n_expression = term_det + term_quad + term_jac

    # The problem asks to output the exact value, which is the symbolic expression.
    # It also asks to "output each number in the final equation", which is achieved
    # by printing the string representation of the symbolic formula.
    print(l_k_n_expression)

if __name__ == '__main__':
    solve_lk_n()