import math

def solve_z_decay():
    """
    Calculates (X_1 * X_2)^-1 based on the Z boson decay formulas.
    """
    # Print the plan and initial values
    print("This script calculates the value of (X_1 * X_2)^-1 for Z boson decay into neutrinos.")
    print("The given constants for massless neutrinos are c_V = 1/2 and c_A = 1/2.\n")
    c_V = 0.5
    c_A = 0.5

    # Step 1: Determine X_2
    print("--- Step 1: Determine X_2 ---")
    print("The decay rate is given by the formula:")
    print("Γ(Z -> f f_bar) = (G_F * m_Z^3 / (12 * sqrt(2) * pi)) * (c_V^2 + c_A^2)")
    
    cv2_plus_ca2 = c_V**2 + c_A**2
    print(f"For neutrinos, c_V^2 + c_A^2 = ({c_V})^2 + ({c_A})^2 = {c_V**2} + {c_A**2} = {cv2_plus_ca2}.")
    
    print("Substituting this value, the decay rate for Z -> nu nu_bar is:")
    print("Γ(Z -> nu nu_bar) = (G_F * m_Z^3 / (12 * sqrt(2) * pi)) * (1/2)")
    print("Γ(Z -> nu nu_bar) = G_F * m_Z^3 / (24 * sqrt(2) * pi)")
    
    print("\nWe are given that Γ(Z -> nu nu_bar) = X_2 * G_F * m_Z^3.")
    print("By comparing the two expressions, we find X_2:")
    print("X_2 = 1 / (24 * sqrt(2) * pi)\n")

    # Step 2: Determine X_1
    print("--- Step 2: Determine X_1 ---")
    print("The decay rate Γ is related to the spin-summed squared amplitude, sum|M|^2, by:")
    print("Γ = (sum|M|^2) / ( (2J+1) * 16 * pi * m_Z )")
    print("For the Z boson, spin J=1, so 2J+1 = 3.")
    print("Γ = (sum|M|^2) / (48 * pi * m_Z)")
    
    print("\nWe assume the problem's |M|^2 corresponds to the spin-summed squared amplitude.")
    print("So, |M|^2 = sum|M|^2 = 48 * pi * m_Z * Γ.")
    print("Substituting our expression for Γ:")
    print("|M|^2 = 48 * pi * m_Z * [G_F * m_Z^3 / (24 * sqrt(2) * pi)]")
    
    print("Simplifying by canceling terms (48*pi / (24*pi)) gives:")
    print("|M|^2 = (2 / sqrt(2)) * G_F * m_Z^4 = sqrt(2) * G_F * m_Z^4")
    
    print("\nWe are given that |M|^2 = X_1 * G_F * m_Z^4.")
    print("By comparing, we find X_1:")
    print("X_1 = sqrt(2)\n")

    # Step 3: Calculate the final result
    print("--- Step 3: Calculate (X_1 * X_2)^-1 ---")
    print("Now we compute the product X_1 * X_2:")
    print("X_1 * X_2 = (sqrt(2)) * (1 / (24 * sqrt(2) * pi))")
    print("The sqrt(2) terms cancel out:")
    print("X_1 * X_2 = 1 / (24 * pi)")
    
    print("\nFinally, we take the inverse:")
    print("(X_1 * X_2)^-1 = 24 * pi")

    # Final numerical calculation
    final_value = 24 * math.pi
    print("\n--- Final Result ---")
    print("The equation for the final answer is:")
    print(f"(X_1 * X_2)^-1 = 24 * pi = {final_value}")

solve_z_decay()
<<<24 * pi>>>