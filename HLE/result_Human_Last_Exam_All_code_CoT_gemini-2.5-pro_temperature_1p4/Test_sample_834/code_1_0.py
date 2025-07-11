import math

def solve_beam_force():
    """
    Calculates the force F required to make the deflection at the beam's end zero.
    """
    # Step 1: Define constants and calculate moments of inertia
    print("--- Step 1: Calculate Cross-Sectional Properties ---")
    
    # Given a = 12^(1/4), so a^4 = 12
    a_4 = 12.0

    # Moments of inertia are calculated using the parallel axis theorem.
    # The formulas for the composite shape are:
    # Iss = (55 * a^4) / 12
    # Izz = (73 * a^4) / 12
    
    Iss = (55 * a_4) / 12
    Izz = (73 * a_4) / 12
    
    print(f"Given a^4 = {a_4}")
    print(f"Moment of inertia about the s-axis, Iss = {Iss}")
    print(f"Moment of inertia about the z-axis, Izz = {Izz}\n")

    # Step 2: Calculate beam parameters L and q0
    print("--- Step 2: Calculate Beam Parameters L and q0 ---")
    
    # L = 30 * Izz / 73
    # q0 = 9 * Iss / 55
    L = (30 * Izz) / 73
    q0 = (9 * Iss) / 55
    
    print(f"Beam length parameter L = (30 * {Izz}) / 73 = {L}")
    print(f"Maximum distributed load q0 = (9 * {Iss}) / 55 = {q0}\n")
    
    # Step 3: Set up deflection equation and solve for F
    print("--- Step 3: Formulate and Solve the Deflection Equation ---")
    
    # The total deflection at x = 3L/2 is the sum of deflection from F and q(x).
    # y_total = y_F + y_q = 0
    # Deflection from F at x=3L/2: y_F = (9 * F * L^3) / (8 * EI)
    # Deflection from q(x) at x=3L/2: y_q = - (13 * q0 * L^4) / (240 * EI)
    
    print("The condition for zero deflection at x = 3L/2 is y_F + y_q = 0.")
    print("This leads to the equation: (9 * F * L^3) / (8 * EI) = (13 * q0 * L^4) / (240 * EI)")
    print("The bending stiffness EI cancels, and we solve for F:")
    print("F = (13 * q0 * L * 8) / (240 * 9)")
    print("F = (13 * q0 * L) / 270\n")

    # Step 4: Substitute numerical values to find F
    print("--- Step 4: Calculate the Final Value of F ---")

    # Final symbolic equation for F
    # F = (13 * q0 * L) / 270
    
    # Substitute values and calculate
    final_F = (13 * q0 * L) / 270
    
    print(f"Substituting the values for q0 and L:")
    print(f"F = (13 * {q0} * {L}) / 270")
    print(f"F = {final_F}")
    
    return final_F

if __name__ == "__main__":
    force_F = solve_beam_force()
    print(f"\nFinal answer for the force F:")
    print(f"<<<{force_F}>>>")