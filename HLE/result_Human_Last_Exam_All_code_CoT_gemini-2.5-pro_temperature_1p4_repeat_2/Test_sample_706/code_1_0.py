from fractions import Fraction

def solve_asymptotic_speed():
    """
    This script calculates the limit of the asymptotic speed of a biased random walk
    on a random ladder graph as the bias c -> infinity.
    """
    # Step 1: Define the probabilities from the problem statement.
    p_v_del = Fraction(1, 2)  # Probability a vertical edge is deleted
    p_h_del = Fraction(1, 3)  # Probability an upper horizontal edge is deleted

    p_v_kept = 1 - p_v_del
    p_h_kept = 1 - p_h_del

    print("Step 1: System Parameters")
    print(f"Probability of vertical edge deletion: {p_v_del}")
    print(f"Probability of upper horizontal edge deletion: {p_h_del}")
    print("-" * 40)

    # Step 2: Analyze the walker's path in the limit c -> infinity.
    # The walker becomes greedy, always choosing the move that maximizes horizontal distance.
    # Priority: Right > Vertical > Left.

    print("Step 2: Limiting Behavior (c -> infinity)")
    print("The walk becomes deterministic, following a greedy path:")
    print("1. Move Right if possible.")
    print("2. Else, move Vertically.")
    print("3. Else, move Left.")
    print("-" * 40)

    # Step 3: Calculate T0, the time to advance one unit from the lower rail.
    # The lower rail is always intact, so a rightward move is always possible.
    # The greedy path is (n, 0) -> (n+1, 0), taking 1 step.
    T0 = Fraction(1)
    print("Step 3: Calculate T0 (Average time from lower rail)")
    print(f"The lower rail is always fully connected. A rightward step is always possible.")
    print(f"It takes 1 step to advance 1 unit. So, T0 = {T0}.")
    print("-" * 40)

    # Step 4: Calculate T1, the expected time to advance one unit from the upper rail.
    # We solve the recursive equation for T1:
    # T1 = p_h_kept * 1 + p_h_del * [p_v_kept * 2 + p_v_del * (1 + 2*T1)]
    # This can be written as T1 = A + B * T1, which gives T1 = A / (1 - B).
    constant_term = p_h_kept + p_h_del * (p_v_kept * 2 + p_v_del)
    T1_coeff = p_h_del * p_v_del * 2
    T1 = constant_term / (1 - T1_coeff)

    print("Step 4: Calculate T1 (Expected time from upper rail)")
    print("We solve a recursive equation for T1 based on the greedy path:")
    print("T1 = P(can_go_right) * 1_step + P(stuck) * E[time_if_stuck]")
    print(f"   = ({p_h_kept})*1 + ({p_h_del})*[({p_v_kept})*2 + ({p_v_del})*(1 + 2*T1)]")
    print(f"Solving this equation for T1 yields T1 = {T1}.")
    print("-" * 40)

    # Step 5: Determine the stationary distribution (pi_0, pi_1).
    # As c -> infinity, transitions L->U are suppressed (~exp(-c)),
    # while transitions U->L are not suppressed.
    # The lower rail acts as a trap, so pi_0 -> 1 and pi_1 -> 0.
    pi_0 = 1
    pi_1 = 0

    print("Step 5: Determine the Stationary Distribution")
    print("In the limit, the lower rail acts as a trap.")
    print("Transitions from lower to upper rail have vanishing probability.")
    print("The stationary distribution for the walker's vertical position is:")
    print(f"pi_0 (on lower rail) = {pi_0}")
    print(f"pi_1 (on upper rail) = {pi_1}")
    print("-" * 40)
    
    # Step 6: Calculate the final asymptotic speed.
    # v = 1 / T_eff = 1 / (pi_0 * T0 + pi_1 * T1)
    T_eff = pi_0 * T0 + pi_1 * T1
    speed = 1 / T_eff

    print("Step 6: Calculate the Final Asymptotic Speed v")
    print("The speed v = 1 / T_effective, where T_effective = pi_0 * T0 + pi_1 * T1.")
    # The problem asks to output the numbers in the final equation.
    final_equation = f"v = 1 / ({pi_0} * {T0} + {pi_1} * {T1})"
    print(f"The final equation is: {final_equation}")
    print(f"T_effective = {pi_0} * {T0} + {pi_1} * {T1} = {T_eff}")
    print(f"v = 1 / {T_eff} = {speed}")
    print("=" * 40)

    # Print the final numerical answer.
    print(f"The final answer is {float(speed)}")

solve_asymptotic_speed()