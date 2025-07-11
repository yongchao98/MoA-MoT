def solve_random_walk():
    """
    Calculates the asymptotic speed of a biased random walk on a random ladder graph
    in the limit of infinite bias (c -> infinity).
    """

    # Step 1: Define probabilities of edge existence from the problem statement.
    # Probability of an upper horizontal edge existing.
    p_h = 1 - 1/3
    # Probability of a vertical edge existing.
    p_v = 1 - 1/2

    print("### Step-by-step Calculation ###")
    print("\nStep 1: Edge Probabilities")
    print(f"The probability of an upper horizontal edge existing is p_h = 1 - 1/3 = {p_h:.4f}")
    print(f"The probability of a vertical edge existing is p_v = 1 - 1/2 = {p_v:.4f}")

    # Step 2: Calculate the limiting speed on the lower level (v0).
    # On level 0, rightward edges always exist. As c->inf, the walker always moves right.
    # Each step gives a displacement of +1.
    v0_limit = 1.0
    print("\nStep 2: Speed on Lower Level (v0)")
    print("On the lower level, rightward edges always exist. The extreme bias forces the walker to always move right.")
    print(f"Therefore, the speed v0 approaches {v0_limit:.4f}")

    # Step 3: Calculate the limiting speed on the upper level (v1).
    # This is the expected displacement per step, averaged over the graph randomness.
    # Displacement is +1 if right edge exists (prob p_h).
    # Displacement is 0 if right edge is missing and vertical exists (prob (1-p_h)*p_v).
    # Displacement is -1 if right and vertical are missing (prob (1-p_h)*(1-p_v)).
    v1_limit = (p_h * 1) + (1 - p_h) * (p_v * 0 + (1 - p_v) * -1)
    print("\nStep 3: Speed on Upper Level (v1)")
    print("On the upper level, the speed is the expected displacement per step:")
    print(f"v1 = (p_h * 1) + (1-p_h) * (p_v * 0 + (1-p_v) * -1)")
    print(f"v1 = ({p_h:.4f} * 1) + ((1-{p_h:.4f}) * ({p_v:.4f} * 0 + (1-{p_v:.4f}) * -1)) = {v1_limit:.4f}")

    # Step 4: Calculate the limiting stationary distribution (pi_0, pi_1).
    # The ratio pi_1/pi_0 = R_0_to_1 / R_1_to_0.
    # R_0_to_1 is the rate of moving up, which is proportional to exp(-c) -> 0.
    # R_1_to_0 is the rate of moving down, which is (1-p_h)*p_v.
    R_0_to_1_limit = 0.0
    R_1_to_0_limit = (1 - p_h) * p_v
    # As R_0_to_1 -> 0, the ratio pi_1/pi_0 -> 0.
    # Since pi_0 + pi_1 = 1, this means pi_0 -> 1 and pi_1 -> 0.
    pi_0_limit = 1.0
    pi_1_limit = 0.0
    print("\nStep 4: Stationary Distribution (pi_0, pi_1)")
    print("The transition rate from level 0 to 1 (R_0_to_1) approaches 0 as c->inf.")
    print("The transition rate from level 1 to 0 (R_1_to_0) approaches a constant.")
    print(f"R_1_to_0 = (1-p_h)*p_v = (1-{p_h:.4f})*{p_v:.4f} = {R_1_to_0_limit:.4f}")
    print("This means the walker spends almost all its time on the lower level.")
    print(f"Therefore, pi_0 approaches {pi_0_limit:.4f} and pi_1 approaches {pi_1_limit:.4f}")

    # Step 5: Calculate the final asymptotic speed v.
    # v = pi_0 * v0 + pi_1 * v1
    v_limit = pi_0_limit * v0_limit + pi_1_limit * v1_limit
    print("\nStep 5: Final Asymptotic Speed (v)")
    print("The overall speed is the weighted average: v = pi_0 * v0 + pi_1 * v1")
    print(f"v = ({pi_0_limit:.4f} * {v0_limit:.4f}) + ({pi_1_limit:.4f} * {v1_limit:.4f})")
    print(f"v = {v_limit:.4f}")

    print("\n" + "="*30)
    print("Final Answer")
    print("="*30)
    print(f"The limit of the asymptotic speed v(c) as c approaches infinity is {v_limit}.")

solve_random_walk()