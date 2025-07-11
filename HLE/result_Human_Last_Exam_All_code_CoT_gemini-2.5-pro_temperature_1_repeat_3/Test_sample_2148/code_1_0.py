import math

def solve_decay_problem():
    """
    Calculates (X1*X2)^-1 based on the provided physics equations.
    """
    # The problem provides the following formulas:
    # 1. Gamma(Z -> f f_bar) = (G_F * m_Z^3) / (12 * sqrt(2) * pi) * (c_V^2 + c_A^2)
    # 2. |M|^2_{m_f=0} = X1 * G_F * m_Z^4
    # 3. Gamma(Z -> nu nu_bar) = X2 * G_F * m_Z^3
    # For neutrinos, c_V = c_A = 1/2.

    # Step 1: Calculate X2.
    # From formula 1 and 3, by comparing the coefficients of G_F * m_Z^3:
    # X2 = (1 / (12 * sqrt(2) * pi)) * (c_V^2 + c_A^2)
    c_V = 0.5
    c_A = 0.5
    # The term (c_V^2 + c_A^2) becomes (0.25 + 0.25) = 0.5
    # So, X2 = 1 / (24 * sqrt(2) * pi)
    X2 = 1 / (24 * math.sqrt(2) * math.pi)

    # Step 2: Calculate X1.
    # The decay rate Gamma is related to the spin-averaged squared amplitude <|M|^2> by:
    # Gamma = <|M|^2> / (16 * pi * m_Z)
    # Assuming |M|^2 in the problem means <|M|^2>, we have:
    # Gamma = (X1 * G_F * m_Z^4) / (16 * pi * m_Z) = (X1 * G_F * m_Z^3) / (16 * pi)
    # Now, we equate this with the expression for Gamma for neutrinos (X2 * G_F * m_Z^3):
    # (X1 * G_F * m_Z^3) / (16 * pi) = X2 * G_F * m_Z^3
    # X1 / (16 * pi) = X2
    # X1 = 16 * pi * X2
    # Substituting the expression for X2:
    # X1 = 16 * pi * (1 / (24 * sqrt(2) * pi))
    # X1 = 16 / (24 * sqrt(2)) = 2 / (3 * sqrt(2)) = sqrt(2) / 3
    X1 = math.sqrt(2) / 3

    # Step 3: Calculate the final result (X1 * X2)^-1
    # product = (sqrt(2) / 3) * (1 / (24 * sqrt(2) * pi)) = 1 / (72 * pi)
    # inverse = 72 * pi
    result = 72 * math.pi

    # Print the numbers for the final equation: (X1 * X2)^-1 = Result
    print("The final equation is of the form: (X1 * X2)^-1 = Result")
    print("The values for each number in the equation are:")
    print(f"X1 = {X1}")
    print(f"X2 = {X2}")
    print(f"Result = {result}")

solve_decay_problem()