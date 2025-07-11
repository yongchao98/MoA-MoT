import math

def solve_random_walk():
    """
    Calculates the limit of the asymptotic speed v(c) as c approaches infinity
    for the described random walk.
    """

    # Let v be the final asymptotic speed.
    # The speed v is the limit as c -> infinity of v(c).
    # The speed v(c) can be written as a weighted average of the speeds on each level:
    # v(c) = pi_0(c) * v_0(c) + pi_1(c) * v_1(c)
    # where:
    # pi_0(c), pi_1(c) are the stationary probabilities of being on the lower {0} and upper {1} levels.
    # v_0(c), v_1(c) are the average speeds on the lower and upper levels for a given c.

    # We need to find the limit of this expression as c -> infinity.

    # Step 1: Determine the limit of the stationary distribution (pi_0, pi_1).

    # The transition from the lower level (0) to the upper level (1) is proportional to e^{c*(x_1_target - x_1_source)} = e^{c*0} = 1.
    # This competes with a move to the right, which is proportional to e^c.
    # The probability of moving up P(0->1) is approximately 1 / e^c, which tends to 0.
    lim_prob_0_to_1 = 0

    # The transition from the upper level (1) to the lower level (0) happens when the walker
    # is blocked horizontally. The probability of an upper horizontal edge being deleted is 1/3.
    # When blocked, the walker prefers moving down (weight 1) to moving left (weight e^-c).
    # Thus, there's a persistent, non-zero probability of moving from the upper to the lower level.
    # Because P(0->1) -> 0 and P(1->0) > 0, the lower level acts as an absorbing state.
    # In the long run, the walker will be on the lower level with probability 1.
    lim_pi_0 = 1
    lim_pi_1 = 1 - lim_pi_0

    # Step 2: Determine the limit of the speed on the lower level, v_0(c).

    # On the lower level, horizontal edges always exist. At any node (n,0), the walker
    # can move to (n+1,0).
    # The probability of moving right is e^c / (e^c + e^-c + ...), which approaches 1 as c -> infinity.
    # A step to the right has a horizontal displacement of +1 and takes 1 time unit.
    # Thus, the speed on the lower level approaches 1.
    displacement_per_step_on_lower = 1
    time_per_step_on_lower = 1
    lim_v_0 = displacement_per_step_on_lower / time_per_step_on_lower

    # Step 3: Determine the limit of the speed on the upper level, v_1(c).

    # The speed on the upper level is bounded (it cannot be more than 1). Its exact
    # value is not needed because it will be multiplied by lim_pi_1, which is 0.
    lim_v_1_desc = "(a bounded value)"


    # Step 4: Calculate the final limit of the speed v.
    # v = lim_{c->inf} v(c) = (lim pi_0) * (lim v_0) + (lim pi_1) * (lim v_1)
    final_v = lim_pi_0 * lim_v_0 + lim_pi_1 * 0 # Multiplying by 0 as lim_v_1 is bounded

    print("The asymptotic speed is given by the formula: v = (lim pi_0) * (lim v_0) + (lim pi_1) * (lim v_1)")
    print("We evaluate each term in the limit c -> infinity:\n")
    print(f"1. The probability of being on the lower level, lim pi_0, is {lim_pi_0}.")
    print(f"2. The probability of being on the upper level, lim pi_1, is {lim_pi_1}.")
    print(f"3. The speed on the lower level, lim v_0, is {lim_v_0}.")
    print(f"4. The speed on the upper level, lim v_1, is {lim_v_1_desc}.")
    print("\nSubstituting these into the formula:\n")
    print(f"v = ({lim_pi_0}) * ({lim_v_0}) + ({lim_pi_1}) * {lim_v_1_desc}")
    print(f"v = {lim_pi_0 * lim_v_0} + 0")
    print(f"v = {final_v}")


solve_random_walk()
<<<1>>>