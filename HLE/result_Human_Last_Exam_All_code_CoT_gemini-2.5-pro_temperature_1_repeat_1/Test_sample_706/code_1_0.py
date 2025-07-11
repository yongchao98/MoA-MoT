import sys

def solve():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    
    The final speed is given by the formula:
    v_inf = pi0_inf * v0_inf + pi1_inf * v1_inf
    
    where:
    - pi0_inf, pi1_inf are the stationary probabilities of being on the lower/upper rail.
    - v0_inf, v1_inf are the asymptotic speeds on the lower/upper rail.
    """

    # --- Step 1: Define edge probabilities ---
    # Probability an upper horizontal edge exists
    p_upper_h = 2/3
    # Probability a vertical edge exists
    p_vertical = 1/2


    # --- Step 2: Determine the limiting stationary distribution ---
    # As c -> infinity, the walker is strongly incentivized to move to the lower rail,
    # which provides an unobstructed path to the right. The probability of jumping
    # up from the lower rail becomes negligible compared to moving right.
    # Therefore, the stationary distribution concentrates entirely on the lower rail.
    pi0_inf = 1.0
    pi1_inf = 0.0


    # --- Step 3: Calculate the limiting speed on the lower rail (v0_inf) ---
    # The lower horizontal edges always exist. As c -> infinity, the walker at (n,0)
    # will always jump to (n+1,0). This is a displacement of +1 in 1 time step.
    v0_inf = 1.0


    # --- Step 4: Calculate the limiting speed on the upper rail (v1_inf) ---
    # This speed is the expected displacement per step, averaged over all possible
    # edge configurations around a vertex on the upper rail, conditional on the
    # vertex being part of the infinite component.
    
    p_R1 = p_upper_h  # Probability of right edge from (n,1)
    p_L1 = p_upper_h  # Probability of left edge from (n,1)
    p_v = p_vertical   # Probability of vertical edge from (n,1)

    # A vertex (n,1) is locally isolated if its right, left, and down edges are all missing.
    # Such a vertex cannot be in the infinite component.
    p_isolated = (1 - p_R1) * (1 - p_L1) * (1 - p_v)
    p_in_component = 1 - p_isolated

    # We calculate the expected displacement, X, by considering the greedy choices of the walker.
    # The final speed v1_inf is E[X | in_component] = E[X * I_in_component] / P(in_component).
    # Since X is undefined for an isolated vertex, its contribution to E[X * I_in_component] is 0.

    # Case 1: Right edge exists. The walker moves right. Displacement = +1.
    # This happens with probability p_R1.
    exp_disp_num = p_R1 * 1
    
    # Case 2: No right edge, but a vertical edge exists. Walker moves down. Displacement = 0.
    # This happens with probability (1 - p_R1) * p_v.
    exp_disp_num += (1 - p_R1) * p_v * 0

    # Case 3: No right or vertical edge, but left edge exists. Walker is forced left. Displacement = -1.
    # This happens with probability (1 - p_R1) * (1 - p_v) * p_L1.
    exp_disp_num += (1 - p_R1) * (1 - p_v) * p_L1 * (-1)
    
    # The speed v1_inf is the conditional expectation.
    v1_inf = exp_disp_num / p_in_component


    # --- Step 5: Combine the results for the final asymptotic speed ---
    v_inf = pi0_inf * v0_inf + pi1_inf * v1_inf
    
    print("The final asymptotic speed is calculated using the equation:")
    print("v_inf = pi0_inf * v0_inf + pi1_inf * v1_inf")
    print("\nIn the limit c -> infinity, the values of these components are:")
    print(f"pi0_inf (stationary prob on lower rail) = {pi0_inf}")
    print(f"v0_inf (speed on lower rail) = {v0_inf}")
    print(f"pi1_inf (stationary prob on upper rail) = {pi1_inf}")
    print(f"v1_inf (speed on upper rail) = {v1_inf:.4f} (which is {exp_disp_num:.4f}/{p_in_component:.4f})")
    
    print("\nPlugging in the numbers:")
    print(f"v_inf = {pi0_inf} * {v0_inf} + {pi1_inf} * {v1_inf:.4f}")
    print(f"v_inf = {v_inf}")

solve()