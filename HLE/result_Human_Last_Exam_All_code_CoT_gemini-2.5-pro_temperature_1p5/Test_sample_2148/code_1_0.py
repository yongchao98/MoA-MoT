import math

def solve_z_decay_problem():
    """
    Calculates (X1 * X2)^-1 based on the Z boson decay formulas.
    """
    # Step 1: Define given values for massless neutrinos
    c_V = 0.5
    c_A = 0.5

    # Calculate the sum of the squared coupling constants
    couplings_sq_sum = c_V**2 + c_A**2

    # Step 2: Determine X2
    # From the problem statement, we compare the two expressions for the decay rate Gamma:
    # Gamma = (G_F * m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    # Gamma = X2 * G_F * m_Z^3
    # Equating them gives: X2 = (c_V^2 + c_A^2) / (12*sqrt(2)*pi)
    X2 = couplings_sq_sum / (12 * math.sqrt(2) * math.pi)

    # Step 3: Determine X1
    # We relate the decay rate Gamma to the squared amplitude |M|^2 using:
    # Gamma = |M|^2 / (48 * pi * m_Z)
    # By equating the two expressions for Gamma, we solve for |M|^2:
    # |M|^2 / (48*pi*m_Z) = (G_F*m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    # |M|^2 = (48*pi*m_Z) * (G_F*m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    # |M|^2 = (2*sqrt(2)) * G_F * m_Z^4 * (c_V^2 + c_A^2)
    # We are given |M|^2 = X1 * G_F * m_Z^4, so by comparing:
    # X1 = 2*sqrt(2) * (c_V^2 + c_A^2)
    X1 = 2 * math.sqrt(2) * couplings_sq_sum

    # Step 4: Calculate the final result
    product = X1 * X2
    final_result = 1 / product
    
    # Step 5: Print the results for each number in the equation
    print(f"Given c_V = {c_V} and c_A = {c_A}")
    print(f"The sum of squared couplings (c_V^2 + c_A^2) is: {couplings_sq_sum}")
    print("-" * 30)
    print(f"The calculated value for X1 is: {X1}")
    print(f"The calculated value for X2 is: {X2}")
    print("-" * 30)
    print(f"The product X1 * X2 is: {product}")
    print(f"The final result for (X1 * X2)^-1 is: {final_result}")

solve_z_decay_problem()
<<<75.39822368615503>>>