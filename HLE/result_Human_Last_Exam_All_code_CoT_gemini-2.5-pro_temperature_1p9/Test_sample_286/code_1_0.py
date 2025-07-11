import math

def solve_for_M():
    """
    This function solves for the smallest positive integer M based on the derived inequality.
    """
    # Let N = 100000.
    # The problem asks for the smallest positive integer M such that for any sequence
    # a_i in [0,1] (for i=1..N), there exists a sequence x_i in [-C, C]
    # (for i=0..N) with C = 1 - 10**-M, satisfying |x_{i-1} - x_i| = a_i.

    # Step 1: Reformulate the condition on x_i.
    # The condition |x_{i-1} - x_i| = a_i means x_i = x_{i-1} +/- a_i.
    # Let's write x_k = x_0 + sum_{i=1 to k} s_i * a_i, where s_i is -1 or 1.
    # Let S_k = sum_{i=1 to k} s_i * a_i. Then x_k = x_0 + S_k.

    # Step 2: Find the condition for the existence of x_0.
    # We need -C <= x_0 + S_k <= C for all k=0..N.
    # This requires that there is an x_0 satisfying max_k(-C - S_k) <= x_0 <= min_k(C - S_k).
    # Such an x_0 exists if and only if max_k(-C - S_k) <= min_k(C - S_k),
    # which simplifies to max_k(S_k) - min_k(S_k) <= 2*C.

    # Step 3: Apply the logic for all `a_i` and some `s_i`, `x_0`.
    # For ANY sequence `a_i`, there must EXIST a choice of signs `s_i`
    # and a starting point `x_0` satisfying the condition.
    # This implies that for any `a_i`, we can choose `s_i` such that the range of partial sums is bounded:
    # min_{s_i} (max_k(S_k) - min_k(S_k)) <= 2*C

    # Step 4: Find the worst-case range.
    # This must hold for the worst-case sequence `a_i`. Let W be this worst-case minimal range.
    # W = max_{a_i} [min_{s_i} (max_k(S_k) - min_k(S_k))]
    # The condition becomes W <= 2*C.
    # A known result in discrepancy theory is that W = 1.
    # A lower bound W >= 1 can be seen by considering a_i = 1 for all i.
    # For a_i=1, the minimal range is achieved by alternating signs, giving partial sums 1,0,1,0,...
    # The range (max - min) for these sums is 1 - 0 = 1.
    # So, we must have 1 <= 2*C.

    # Step 5: Solve for M.
    # We substitute C = 1 - 10**(-M) into the inequality.
    print("The final inequality to solve is: 1 <= 2 * (1 - 10**-M)")
    # 1 / 2 <= 1 - 10**-M
    print("This simplifies to: 0.5 <= 1 - 10**-M")
    # 10**-M <= 1 - 0.5
    print("Which becomes: 10**-M <= 0.5")
    # Taking log10 on both sides: -M <= log10(0.5)
    # -M <= -log10(2)
    # M >= log10(2)
    log10_of_2 = math.log10(2)
    print(f"The condition on M is: M >= log10(2)")
    print(f"Numerically, M >= {log10_of_2}")

    # Since M must be the smallest positive integer satisfying this, we take the ceiling.
    M = math.ceil(log10_of_2)
    print(f"The smallest positive integer M is therefore {M}.")

solve_for_M()
