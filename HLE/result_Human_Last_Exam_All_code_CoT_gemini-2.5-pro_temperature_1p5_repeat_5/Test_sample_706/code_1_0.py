import math

def solve():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.

    The asymptotic speed can be expressed as:
    v(c) = pi_0(c) * Delta_0(c) + pi_1(c) * Delta_1(c)
    where pi_i(c) is the stationary probability of being on level i,
    and Delta_i(c) is the expected horizontal displacement in one step from level i.
    We need to find the limit of v(c) as c -> infinity.

    1. Limit of Stationary Probabilities (pi_0, pi_1):
    In steady state, the probability flow from level 0 to 1 must equal the flow from 1 to 0.
    The probability of a jump 0->1 is proportional to exp(-c), while the probability of a jump 1->0
    is significant (it happens when the rightward path on level 1 is blocked).
    This implies that as c -> infinity, the walker is infinitely more likely to be on level 0.
    lim_{c->inf} pi_0(c) = 1
    lim_{c->inf} pi_1(c) = 0
    """
    pi_0_limit = 1
    pi_1_limit = 0

    """
    2. Limit of Expected Displacements (Delta_0, Delta_1):

    From level 0, at (n,0):
    There is always a rightward edge to (n+1,0). As c->inf, the probability of taking this edge
    approaches 1. The displacement is +1.
    lim_{c->inf} Delta_0(c) = 1

    From level 1, at (n,1):
    A rightward edge to (n+1,1) exists with probability p_h_exist = 1 - 1/3 = 2/3.
    If it exists, the walker moves right with displacement +1.
    If it doesn't exist (prob 1/3), the walker cannot move right. The next best options (vertical
    or left) have non-positive displacement. As c->inf, the horizontal displacement in this case is 0.
    So, the expected displacement is P(right edge exists) * 1 + P(right edge missing) * 0.
    lim_{c->inf} Delta_1(c) = (2/3) * 1 + (1/3) * 0 = 2/3
    """
    p_upper_horizontal_exists = 2/3
    delta_0_limit = 1.0
    delta_1_limit = p_upper_horizontal_exists * 1.0 + (1 - p_upper_horizontal_exists) * 0.0

    # The overall speed limit is the sum of these limiting components.
    # v = (lim pi_0) * (lim Delta_0) + (lim pi_1) * (lim Delta_1)
    v_limit = pi_0_limit * delta_0_limit + pi_1_limit * delta_1_limit

    print("The limit of the asymptotic speed v(c) is calculated as:")
    print("v = (lim pi_0) * (lim Delta_0) + (lim pi_1) * (lim Delta_1)")
    print(f"where:")
    print(f"lim pi_0 = {pi_0_limit} (stationary probability of being on the lower level)")
    print(f"lim pi_1 = {pi_1_limit} (stationary probability of being on the upper level)")
    print(f"lim Delta_0 = {delta_0_limit} (avg. displacement from the lower level)")
    print(f"lim Delta_1 = {delta_1_limit:.4f} (avg. displacement from the upper level)")
    print("\nFinal Equation:")
    print(f"v = ({pi_0_limit}) * ({delta_0_limit}) + ({pi_1_limit}) * ({delta_1_limit:.4f})")
    print(f"v = {v_limit}")

solve()
<<<1>>>