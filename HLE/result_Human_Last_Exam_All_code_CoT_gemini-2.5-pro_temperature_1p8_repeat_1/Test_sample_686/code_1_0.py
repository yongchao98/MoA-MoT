import sympy as sp

def solve_ising_susceptibility():
    """
    This function prints the derived formula for the magnetic susceptibility chi.
    The derivation follows the plan outlined above.
    """

    # Define the symbols used in the equation
    chi, N, c, A = sp.symbols('chi N c A')
    
    # The derivation steps outlined in the plan lead to the following relation:
    # 1. C_l = (1 - m_0^2) * A^l
    # 2. chi = beta * c * (1 - m_0^2) * Sum_{l=1 to inf} ( (c-1) * A )^(l-1) * A
    # 3. Sum = A / (1 - (c-1)*A)
    # 4. chi = beta * c * (1 - m_0^2) * A / (1 - (c-1)*A)
    # 5. Using N = beta * c * (1 - m_0^2) / (c-1), we have beta * c * (1 - m_0^2) = N * (c-1)
    # 6. chi = N * (c-1) * A / (1 - (c-1)*A)

    # Final equation
    equation = sp.Eq(chi, N * ((c - 1) * A) / (1 - (c - 1) * A))
    
    # The task asks to output each number in the final equation.
    # We can rewrite the equation as: chi * (1 - (c-1)*A) - N*(c-1)*A = 0
    # The coefficients of the terms are polynomials in the symbols.
    # The numerical integer constants are 1 and -1.
    
    print("The final equation for the magnetic susceptibility is:")
    print(equation)
    
    # Let k = (c-1)*A. The equation is chi = N*k / (1-k).
    # The explicit integer numbers in this form of the equation are 1 and -1.
    print("\nEach number in the final equation:")
    print("Numerator term (c - 1): coefficient is 1")
    print("Denominator term 1: coefficient is 1")
    print("Denominator term -(c - 1): coefficient is -1")

solve_ising_susceptibility()
