import math

# Step 1: Define a function to calculate geometric properties I_zz and I_ss.
# The cross-section is a 3a x 3a square with two a x a squares removed.
# The centers of the cutouts are at (a/2, -a) and (-a/2, a).
# The value 'a' is given as 12^(1/4), so a^4 = 12.
def calculate_geometric_properties():
    """
    Calculates the second moments of area I_zz and I_ss.
    """
    a_fourth = 12.0
    
    # Moment of inertia for a square b x b about its centroidal axis is b^4 / 12.
    
    # For the large square (3a x 3a):
    I_zz_large = (3**4 * a_fourth) / 12
    I_ss_large = (3**4 * a_fourth) / 12
    
    # For the cutouts (a x a):
    # I_c = a^4 / 12
    # Area = a^2
    # a^2 = (12^(1/4))^2 = 12^(1/2)
    # We can work with a_fourth directly.
    I_c_cutout = a_fourth / 12
    
    # Cutout 1 at (s=a/2, z=-a)
    # d_s1 = a/2, d_z1 = -a
    I_zz_cutout1 = I_c_cutout + (a_fourth**0.5) * (a_fourth**0.25 / 2)**2 # A * d_s^2
    I_ss_cutout1 = I_c_cutout + (a_fourth**0.5) * (-a_fourth**0.25)**2 # A * d_z^2
    # A cleaner way using a^4 and the parallel axis theorem:
    # I_zz = I_c + A*d_s^2 = a^4/12 + a^2*(a/2)^2 = a^4/12 + a^4/4 = a^4/3
    I_zz_cutout_total = 2 * (a_fourth / 3)
    # I_ss = I_c + A*d_z^2 = a^4/12 + a^2*a^2 = a^4/12 + a^4 = 13*a^4/12
    I_ss_cutout_total = 2 * (13 * a_fourth / 12)

    # Total I_zz and I_ss by subtraction
    I_zz = I_zz_large - I_zz_cutout_total
    I_ss = I_ss_large - I_ss_cutout_total

    return I_zz, I_ss

# Step 2: Use the derived relation for F and the given data to find the result.
def solve_force():
    """
    Calculates L, q0, and the final force F.
    """
    I_zz, I_ss = calculate_geometric_properties()
    
    print("Step 1: Calculating geometric properties")
    print(f"Given a^4 = 12")
    print(f"I_zz = 81*a^4/12 - 2*(a^4/3) = 73*a^4/12 = 73*12/12 = {I_zz}")
    print(f"I_ss = 81*a^4/12 - 2*(13*a^4/12) = 55*a^4/12 = 55*12/12 = {I_ss}")
    print("-" * 30)

    # Calculate L and q0 from the given data
    L = (30.0 * I_zz) / 73.0
    q0 = (9.0 * I_ss) / 55.0

    print("Step 2: Calculating L and q0")
    print(f"L = (30 * I_zz) / 73 = (30 * {I_zz:.1f}) / 73 = {L:.1f}")
    print(f"q0 = (9 * I_ss) / 55 = (9 * {I_ss:.1f}) / 55 = {q0:.1f}")
    print("-" * 30)
    
    # The deflection at x=3L/2 is set to zero: y_F + y_q = 0
    # y_F = F * (3L/2)^3 / (3*EI) = 9*F*L^3 / (8*EI)
    # y_q = -13*q0*L^4 / (240*EI)
    # Setting them equal and solving for F:
    # 9*F*L^3 / (8*EI) = 13*q0*L^4 / (240*EI)
    # F = (13 * q0 * L / 240) * (8 / 9)
    # F = 13 * q0 * L / 270
    
    print("Step 3: Deriving and calculating the force F")
    print("The force F is found from the zero-deflection condition:")
    print("F = (13 * q0 * L) / 270")
    print("Substituting the values for q0 and L:")
    
    # Print the final equation with numbers
    print(f"F = (13 * {q0:.1f} * {L:.1f}) / 270")
    
    F = (13 * q0 * L) / 270.0
    
    print(f"F = {F:.1f}")
    
    return F

# Execute the solver and print the final answer in the required format
final_force = solve_force()
print(f"<<<{final_force:.0f}>>>")
