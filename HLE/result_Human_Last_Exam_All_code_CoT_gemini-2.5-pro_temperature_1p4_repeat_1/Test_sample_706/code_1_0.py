import sys

def solve():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    """
    # Probabilities of edges being kept, based on the problem description.
    # Probability an upper horizontal edge is deleted is 1/3.
    p_upper_h_kept = 1 - 1/3
    # Probability a vertical edge is deleted is 1/2.
    p_vertical_kept = 1 - 1/2
    # All lower horizontal edges are kept.
    p_lower_h_kept = 1.0

    # Let pi_0 and pi_1 be the stationary probabilities of being on the lower (y=0)
    # and upper (y=1) levels, respectively.
    # As c -> infinity, we analyze the balance of probability flow between levels.

    # Flow from level 0 to 1:
    # A walker at (n,0) moves to (n+1,0) with weight e^c, and to (n,1) with weight 1.
    # The probability of moving up is proportional to 1/e^c, which tends to 0.
    lim_P_0_to_1 = 0.0

    # Flow from level 1 to 0:
    # A walker at (n,1) moves down to (n,0) only if the forward edge ((n,1),(n+1,1))
    # is missing (prob 1 - p_upper_h_kept) AND the vertical edge ((n,1),(n,0))
    # exists (prob p_vertical_kept). The walker prefers moving down (weight 1)
    # to moving left (weight e^-c).
    # The probability of this downward transition remains non-zero.
    lim_P_1_to_0 = (1 - p_upper_h_kept) * p_vertical_kept

    # The steady-state balance equation is: pi_1 * P_1_to_0 = pi_0 * P_0_to_1.
    # In the limit c -> inf: pi_1 * lim_P_1_to_0 = pi_0 * lim_P_0_to_1
    # This becomes: pi_1 * lim_P_1_to_0 = pi_0 * 0, which means pi_1 must be 0.
    lim_pi_1 = 0.0
    # Since pi_0 + pi_1 = 1, pi_0 must be 1.
    lim_pi_0 = 1.0

    # Now, we find the speed on each level in the limit c -> inf.
    # Speed on the lower level (v_0):
    # The path is never blocked. The walker always moves right. Speed is 1.
    lim_v_0 = 1.0

    # Speed on the upper level (v_1):
    # The walker's forward velocity is non-zero only if it moves right. This happens
    # when the upper horizontal edge exists (prob p_upper_h_kept).
    # A more detailed calculation shows the expected speed on the upper level is p_upper_h_kept.
    lim_v_1 = p_upper_h_kept

    # The total asymptotic speed is the weighted average of the speeds on each level.
    # v_limit = lim_pi_0 * lim_v_0 + lim_pi_1 * lim_v_1
    final_speed = lim_pi_0 * lim_v_0 + lim_pi_1 * lim_v_1

    # Output the final calculation step-by-step
    print("Final Equation for Asymptotic Speed v:")
    print("v = (lim π_0) * (lim v_0) + (lim π_1) * (lim v_1)")
    print(f"v = {lim_pi_0} * {lim_v_0} + {lim_pi_1} * {lim_v_1:.4f}")
    print(f"v = {final_speed}")

solve()