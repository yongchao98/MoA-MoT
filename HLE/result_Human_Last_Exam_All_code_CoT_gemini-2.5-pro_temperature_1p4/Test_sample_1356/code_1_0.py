import math

def solve_pendulum_period():
    """
    This function calculates and explains the derivation for the period of motion
    of the described two-disk pendulum system.
    """
    print("Derivation of the period of motion for the two-disk system:")
    print("Let m be the mass and R be the radius of each disk.")
    print("Let g be the acceleration due to gravity.")
    print("Let x be the horizontal position of the top disk's center, and θ be the angle of the rod with the vertical.\n")

    print("Step 1: Define the Kinetic (T) and Potential (V) energies for small angles.")
    print(" - Kinetic energy of the top disk (rolling without slipping): T1 = (3/4)*m*(dx/dt)^2")
    print(" - Kinetic energy of the bottom disk (translation only): T2 ≈ (1/2)*m*(dx/dt + 4*R*dθ/dt)^2")
    print(" - Total Kinetic Energy T = T1 + T2 ≈ (5/4)*m*(dx/dt)^2 + 4*m*R*(dx/dt)*(dθ/dt) + 8*m*R^2*(dθ/dt)^2")
    print(" - Potential Energy V (relative to the top disk's center) ≈ 2*m*g*R*θ^2\n")

    print("Step 2: Find the constraint between x and θ using the Lagrangian.")
    print("The Lagrangian L = T - V does not depend on the position x, so the momentum p_x is conserved.")
    print("p_x = ∂L/∂(dx/dt) = (5/2)*m*(dx/dt) + 4*m*R*(dθ/dt)")
    print("For a system starting from rest, p_x = 0, which gives the constraint:")
    print("dx/dt = -(8/5)*R*(dθ/dt)\n")

    print("Step 3: Create an effective Lagrangian in terms of θ only.")
    print("Substituting the constraint into T and V gives:")
    print(" - Effective Kinetic Energy: T_eff = (24/5)*m*R^2*(dθ/dt)^2")
    print(" - Effective Potential Energy: V_eff = 2*m*g*R*θ^2\n")

    print("Step 4: Form the equation of motion for small oscillations.")
    print("The equation for a simple harmonic oscillator is M_eff * d^2θ/dt^2 + K_eff * θ = 0.")
    print("Comparing T_eff = (1/2)*M_eff*(dθ/dt)^2 and V_eff = (1/2)*K_eff*θ^2, we get:")
    print(" - Effective 'mass' M_eff = (48/5)*m*R^2")
    print(" - Effective 'spring constant' K_eff = 4*m*g*R\n")

    print("Step 5: Calculate the angular frequency (ω) and the period (T).")
    print("The angular frequency squared is ω^2 = K_eff / M_eff.")
    print("ω^2 = (4*m*g*R) / ((48/5)*m*R^2) = (4 * 5) / 48 * (g/R) = 20/48 * (g/R) = 5/12 * (g/R)")
    print("The period T = 2*π / ω.\n")

    print("--- Final Answer ---")
    numerator = 12
    denominator = 5
    print("The final equation for the period T is:")
    print(f"T = 2 * π * sqrt( ({numerator} * R) / ({denominator} * g) )")

solve_pendulum_period()