from fractions import Fraction

def solve_random_walk_speed():
    """
    Calculates the limit of the asymptotic speed v(c) as c approaches infinity for a biased random walk on a random ladder graph.

    The method proceeds in several steps:
    1.  Define the probabilities for edges being kept in the random graph.
    2.  Express the overall speed v(c) as a weighted average of the speeds on the lower (v_0) and upper (v_1) rails, weighted by the stationary probabilities (pi_0, pi_1) of being on each rail.
        v(c) = pi_0(c) * v_0(c) + pi_1(c) * v_1(c)
    3.  Analyze the system in the limit c -> infinity. We find the limits of v_0, v_1, pi_0, and pi_1.
    4.  The speed on a rail is the expected horizontal displacement per step for a walker on that rail.
    5.  The stationary probabilities are found by balancing the probability flow between the two rails in the steady state.
    6.  The final result is obtained by combining these limiting values.
    """

    # --- Step 1: Define edge probabilities from the problem description ---
    # Probability a vertical edge is deleted is 1/2.
    p_v_kept = Fraction(1, 2)
    # Probability an upper horizontal edge is deleted is 1/3.
    p_hu_kept = Fraction(2, 3)

    # --- Step 2: Calculate lim v_0(c), the speed on the lower rail ---
    # On the lower rail, the rightward horizontal edge always exists. As c -> infinity,
    # the jump weight e^c for moving right massively outweighs the weight for any other
    # move. Therefore, the walker almost always moves right. Each step takes 1 unit of time
    # and provides +1 horizontal displacement.
    v0 = Fraction(1)

    # --- Step 3: Calculate lim v_1(c), the speed on the upper rail ---
    # The speed on the upper rail depends on the local edge configuration. We average
    # the limiting speed over all configurations, conditioning on the node not being isolated.
    # Let R, L, V be indicators for Right, Left, Vertical edges existing.
    # Limiting speed contributions: +1 if R=1; -1 if R=0,V=0,L=1; 0 if R=0,V=1.
    p_R1 = p_hu_kept
    p_R0_V1 = (1 - p_hu_kept) * p_v_kept
    p_R0_V0_L1 = (1 - p_hu_kept) * (1 - p_v_kept) * p_hu_kept
    # The probability that a node is isolated from its three neighbors (R=0,L=0,V=0):
    p_isolated = (1 - p_hu_kept) * (1 - p_v_kept) * (1 - p_hu_kept)
    # We condition on the node not being isolated by normalizing.
    p_norm = 1 - p_isolated
    # The average speed is the sum of (speed * probability), normalized.
    v1_numerator = (p_R1 * 1) + (p_R0_V1 * 0) + (p_R0_V0_L1 * -1)
    v1 = v1_numerator / p_norm

    # --- Step 4: Calculate lim pi_0(c) and lim pi_1(c) ---
    # In equilibrium, probability flow is balanced: pi_0 * P(0->1) = pi_1 * P(1->0).
    # We examine the average transition probabilities in the c -> inf limit.
    # P(0->1): Jump up (weight 1) vs jump right (weight e^c). lim E[P(0->1)] = 0.
    lim_p_01 = Fraction(0)
    # P(1->0): If the right edge is missing, jump down (weight 1) competes with
    # jump left (weight e^-c). The jump down is overwhelmingly preferred.
    # lim E[P(1->0)] = P(Right edge missing AND Vertical edge exists)
    lim_p_10 = (1 - p_hu_kept) * p_v_kept
    # The balance equation becomes: pi_0 * 0 = pi_1 * (1/6).
    # Since lim_p_10 is non-zero, this implies lim pi_1 = 0.
    pi0 = Fraction(1)
    pi1 = Fraction(0)

    # --- Step 5: Final Result ---
    v_final = pi0 * v0 + pi1 * v1

    print("Calculation of the asymptotic speed v = lim_{c->inf} v(c)")
    print("---------------------------------------------------------")
    print("The speed is v = pi_0 * v_0 + pi_1 * v_1, where all values are limits as c -> inf.")
    print(f"\n1. Limiting speed on lower rail, v_0:\n   v_0 = {v0}")
    print(f"\n2. Limiting speed on upper rail, v_1:\n   v_1 = {v1_numerator}/{p_norm} = {v1}")
    print(f"\n3. Limiting stationary probabilities, (pi_0, pi_1):")
    print(f"   From the balance equation pi_0 * lim P(0->1) = pi_1 * lim P(1->0)")
    print(f"   pi_0 * {lim_p_01} = pi_1 * {lim_p_10}")
    print(f"   This implies pi_1 = {pi1}, and therefore pi_0 = {pi0}.")
    print("\n4. Final Result Equation:")
    print(f"   v = ({pi0}) * ({v0}) + ({pi1}) * ({v1})")
    print(f"   v = {v_final}")

solve_random_walk_speed()
<<<1>>>