import math

def solve_z_decay():
    """
    Calculates (X1*X2)^-1 based on the Z boson decay equations.
    """
    # Step 1: Define constants for the neutrino case
    c_V = 0.5
    c_A = 0.5
    
    # This combination appears in both X1 and X2
    cv2_plus_ca2 = c_V**2 + c_A**2
    
    print(f"For neutrinos, c_V = {c_V} and c_A = {c_A}.")
    print(f"The term (c_V^2 + c_A^2) is calculated as: {cv2_plus_ca2}\n")

    # Step 2: Calculate X1
    # From the relationship Gamma = |M|^2 / (16*pi*m_Z), we derive
    # |M|^2 = 16*pi*m_Z * Gamma.
    # Substituting Gamma = (G_F*m_Z^3)/(12*sqrt(2)*pi) * (c_V^2+c_A^2), we get:
    # |M|^2 = (4/(3*sqrt(2))) * G_F*m_Z^4 * (c_V^2+c_A^2)
    # Comparing with |M|^2 = X1*G_F*m_Z^4, we find X1.
    X1 = (4 / (3 * math.sqrt(2))) * cv2_plus_ca2
    
    # Step 3: Calculate X2
    # Comparing Gamma = (G_F*m_Z^3)/(12*sqrt(2)*pi) * (c_V^2+c_A^2)
    # with Gamma = X2*G_F*m_Z^3, we find X2.
    X2 = cv2_plus_ca2 / (12 * math.sqrt(2) * math.pi)
    
    # Step 4: Calculate the final result
    product_X1_X2 = X1 * X2
    final_result = 1 / product_X1_X2
    
    # Step 5: Print the final equation with all the numbers
    print("Derived expressions:")
    print(f"X1 = (4 / (3*sqrt(2))) * (c_V^2 + c_A^2) = {X1}")
    print(f"X2 = (c_V^2 + c_A^2) / (12*sqrt(2)*pi) = {X2}\n")

    print("Final Calculation:")
    print(f"The final equation is (X1 * X2)^-1")
    print(f"Substituting the values: ( {X1} * {X2} )^-1")
    print(f"= ( {product_X1_X2} )^-1")
    print(f"= {final_result}")

solve_z_decay()