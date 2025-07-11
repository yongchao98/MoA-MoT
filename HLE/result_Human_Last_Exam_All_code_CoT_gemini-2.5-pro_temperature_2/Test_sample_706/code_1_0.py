import numpy as np

def solve():
    """
    This script calculates the asymptotic speed of the random walk by analyzing
    the system's behavior as the bias constant c approaches infinity.
    """
    
    # Probability that a vertical edge exists
    p_v = 1/2
    # Probability that an upper horizontal edge exists
    p_u = 2/3
    
    print(f"Probabilities are p_v = {p_v} and p_u = {p_u:.4f}\n")

    # We use the limiting values directly, as derived in the explanation.
    # Rate of transitioning from lower (0) to upper (1) rail.
    # R_01(c) = p_v / (e^c + e^-c + 1)
    R_01_lim = 0
    
    # Rate of transitioning from upper (1) to lower (0) rail.
    # R_10(c) -> p_v * (1 - p_u) as c -> inf
    R_10_lim = p_v * (1 - p_u)
    
    print("Step 1: Calculate stationary distribution (pi_0, pi_1)")
    print("From the balance equation pi_0 * R_01 = pi_1 * R_10:")
    print(f"lim R_01(c) = {R_01_lim}")
    print(f"lim R_10(c) = {R_10_lim:.4f}")
    
    # pi_1 = R_01 / (R_10 + R_01)
    # As R_01 -> 0 and R_10 -> constant > 0, pi_1 must go to 0.
    pi_1_lim = 0
    pi_0_lim = 1
    
    print(f"--> Probability on upper rail, pi_1, approaches {pi_1_lim}")
    print(f"--> Probability on lower rail, pi_0, approaches {pi_0_lim}\n")
    
    print("Step 2: Calculate speeds on each level (v_0, v_1)")
    # Speed on lower rail, v_0, approaches 1.
    v_0_lim = 1.0
    # Speed on upper rail, v_1, approaches 4/7.
    # This value is not needed for the final result but is part of the system dynamics.
    v_1_lim = 4/7

    print(f"Speed on lower rail, v_0, approaches {v_0_lim}")
    print(f"Speed on upper rail, v_1, approaches {v_1_lim:.4f}\n")
    
    print("Step 3: Combine to find the final asymptotic speed v")
    print("The formula is: v = pi_0 * v_0 + pi_1 * v_1")
    
    final_v = pi_0_lim * v_0_lim + pi_1_lim * v_1_lim

    print("\nSubstituting the limiting values:")
    print(f"v = {pi_0_lim} * {v_0_lim} + {pi_1_lim} * {v_1_lim:.4f}")
    print(f"v = {final_v}")
    
solve()