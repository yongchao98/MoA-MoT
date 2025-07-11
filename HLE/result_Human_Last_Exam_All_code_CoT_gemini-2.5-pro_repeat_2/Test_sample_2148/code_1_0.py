import math

def solve_z_decay():
    """
    Calculates the value of (X1 * X2)^-1 based on the provided physics formulas.
    """
    # Step 1: Define the given coupling constants for massless neutrinos.
    c_V = 0.5
    c_A = 0.5
    print(f"The given coupling constants are: c_V = {c_V}, c_A = {c_A}")

    # Step 2: Calculate the term (c_V^2 + c_A^2), which appears in the formulas.
    C_sum_sq = c_V**2 + c_A**2
    print(f"The term (c_V^2 + c_A^2) evaluates to: {C_sum_sq}")

    # Step 3: Derive the analytical expression for (X1 * X2)^-1.
    # From the problem statement, we are given:
    # 1) |M|^2 = X1 * G_F * m_Z^4
    # 2) Gamma = X2 * G_F * m_Z^3
    # 3) Gamma = (G_F * m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    # The general relation between the spin-summed |M|^2 and Gamma for a 2-body decay is:
    # Gamma = |M|^2 / (48 * pi * m_Z)
    #
    # From this, we derive expressions for X1 and X2:
    # |M|^2 = 48 * pi * m_Z * Gamma
    # X1 * G_F * m_Z^4 = 48 * pi * m_Z * (G_F * m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    # X1 = (48*pi / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2) = 2*sqrt(2) * (c_V^2 + c_A^2)
    #
    # X2 = (1 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
    #
    # Now, we find the product X1 * X2:
    # X1 * X2 = [2*sqrt(2)*(c_V^2+c_A^2)] * [(c_V^2+c_A^2)/(12*sqrt(2)*pi)]
    # X1 * X2 = (1 / (6*pi)) * (c_V^2 + c_A^2)^2
    #
    # The inverse is (X1 * X2)^-1 = 6*pi / (c_V^2 + c_A^2)^2
    # Substituting C_sum_sq = 0.5, we get:
    # (X1 * X2)^-1 = 6*pi / (0.5)^2 = 6*pi / 0.25 = 24*pi.

    # Step 4: Calculate the final numerical result.
    coefficient = 24
    final_result = coefficient * math.pi
    
    print("\nFrom the physics relations, we find the final equation for the result:")
    print(f"(X1 * X2)^-1 = {coefficient} * pi")
    print("\nEvaluating this expression numerically:")
    print(f"{coefficient} * {math.pi} = {final_result}")

solve_z_decay()