import math

def solve():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.

    The method involves analyzing the stationary distribution of the walker
    between the upper and lower rails of the ladder graph in the c -> infinity limit.
    """

    # Probability of an upper horizontal edge being present
    p_u = 2/3
    # Probability of a vertical edge being present
    p_v = 1/2

    print("Step 1: Analyze the stationary distribution as c -> infinity.")
    print("Let pi_0(c) and pi_1(c) be the stationary probabilities of being on the lower and upper rail.")
    print("We balance the probability flow between the rails.")

    # r_01 is the rate of moving from lower to upper rail.
    # This move competes with a rightward move of weight e^c, so r_01 is proportional to e^(-c).
    # r_10 is the rate of moving from upper to lower rail.
    # This happens when a gap is encountered on the upper rail.
    # The rate is proportional to the probability of a gap (1 - p_u).
    
    print(f"The ratio pi_1(c) / pi_0(c) is proportional to r_01 / r_10 ~ e^(-c) / (1 - p_u).")
    print("As c -> infinity, e^(-c) -> 0.")
    
    # As c -> infinity, the ratio pi_1/pi_0 goes to 0.
    pi_1_limit = 0
    # Since pi_0 + pi_1 = 1, pi_0 must go to 1.
    pi_0_limit = 1
    
    print(f"\nTherefore, in the limit c -> infinity:")
    print(f"pi_0 = {pi_0_limit}")
    print(f"pi_1 = {pi_1_limit}")
    
    print("\nStep 2: Calculate the asymptotic speed v.")
    print("The speed is the expected horizontal displacement per step in the stationary state.")
    print("v = pi_0 * E[dx | on lower rail] + pi_1 * E[dx | on upper rail]")
    
    # Expected displacement on the lower rail (always moves right)
    E_dx_0 = 1
    
    # Expected displacement on the upper rail (moves right with prob p_u, down with prob 1-p_u)
    E_dx_1 = p_u * 1 + (1 - p_u) * 0
    
    print(f"\nE[dx | on lower rail] = {E_dx_0}")
    print(f"E[dx | on upper rail] = p_u * 1 + (1 - p_u) * 0 = {p_u:.4f}")
    
    # Final speed calculation
    v_limit = pi_0_limit * E_dx_0 + pi_1_limit * E_dx_1

    print("\nFinal equation for the speed v in the limit:")
    # print each number in the final equation
    print(f"v = ({pi_0_limit}) * ({E_dx_0}) + ({pi_1_limit}) * ({E_dx_1:.4f})")
    print(f"v = {v_limit}")

solve()