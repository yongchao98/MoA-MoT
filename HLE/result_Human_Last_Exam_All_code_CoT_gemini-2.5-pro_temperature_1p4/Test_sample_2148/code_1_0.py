import math

def solve_z_decay():
    """
    Calculates (X_1 * X_2)^(-1) based on the Z boson decay formulas.
    """
    # Given values for massless neutrinos
    c_V = 0.5
    c_A = 0.5

    # Step 1: Calculate the term (c_V^2 + c_A^2)
    cv2_plus_ca2 = c_V**2 + c_A**2
    # cv2_plus_ca2 = (0.5)^2 + (0.5)^2 = 0.25 + 0.25 = 0.5

    # Step 2: Calculate X_2
    # From Γ = (G_F * m_Z^3) / (12 * sqrt(2) * π) * (c_V^2 + c_A^2)
    # and Γ = X_2 * G_F * m_Z^3
    # We get X_2 = (c_V^2 + c_A^2) / (12 * sqrt(2) * π)
    X_2 = cv2_plus_ca2 / (12 * math.sqrt(2) * math.pi)

    # Step 3: Calculate X_1
    # We relate Γ and |M|^2 using the phase space formula: Γ = |M|^2 / (48 * π * m_Z)
    # So, |M|^2 = 48 * π * m_Z * Γ
    # Substituting the definitions:
    # X_1 * G_F * m_Z^4 = 48 * π * m_Z * (X_2 * G_F * m_Z^3)
    # X_1 = 48 * π * X_2
    # Alternatively, from |M|^2 = 2 * sqrt(2) * G_F * m_Z^4 * (c_V^2 + c_A^2)
    # We get X_1 = 2 * sqrt(2) * (c_V^2 + c_A^2)
    X_1 = 2 * math.sqrt(2) * cv2_plus_ca2

    # Step 4: Calculate the final result
    # (X_1 * X_2)^(-1)
    result = 1 / (X_1 * X_2)
    
    # Print the equation with the calculated numbers
    print(f"For c_V = {c_V} and c_A = {c_A}:")
    print(f"The value of X_1 is: {X_1}")
    print(f"The value of X_2 is: {X_2}")
    print(f"The equation is ({X_1}) * ({X_2}) = {X_1 * X_2}")
    print(f"The final result for (X_1 * X_2)^-1 is: {result}")

solve_z_decay()