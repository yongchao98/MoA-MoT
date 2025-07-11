from fractions import Fraction

def solve_committee_problem():
    """
    Solves the committee election problem by calculating s(N,W1) and s(N,W2)
    and their ratio.
    """
    # 1. Define candidates, voters, and approval sets
    C_common = set(range(1, 9))
    P = {i: {8 + 2 * i - 1, 8 + 2 * i} for i in range(1, 9)}
    Q_9 = set(range(25, 29))
    Q_10 = set(range(29, 33))

    A = {i: C_common.union(P[i]) for i in range(1, 9)}
    A[9] = Q_9
    A[10] = Q_10

    N = set(range(1, 9))
    k = 20
    n = 10

    def calculate_satisfaction(group_N, committee_W, approvals_A):
        satisfaction = 0
        for i in group_N:
            satisfaction += len(approvals_A[i].intersection(committee_W))
        return satisfaction

    # 2. Determine s(N, W1) for the core committee
    # A committee W is in the core if for any c_in not in W and c_out in W,
    # it is not the case that the set of approvers of c_out is a strict
    # subset of the approvers of c_in.
    # Voters approving any candidate in P[i] is {i}.
    # Voters approving any candidate in C_common is {1,...,8}.
    # Since {i} is a strict subset of {1,...,8}, any candidate from C_common
    # has priority over any candidate from P[i].
    # Therefore, a core committee must include all of C_common before any of P_i.
    # To minimize s(N,W), we build W1 by:
    # - Including all 8 candidates from C_common.
    # - Including candidates not approved by N, i.e., from Q_9 and Q_10 (8 candidates).
    # - Filling the remaining 4 spots from candidates in P_i sets.
    
    W1_C_common = C_common
    W1_Q = Q_9.union(Q_10)
    # To fill the last 4 spots, we can take P[1] and P[2]
    W1_P_part = P[1].union(P[2])
    W1 = W1_C_common.union(W1_Q).union(W1_P_part)

    s_N_W1 = calculate_satisfaction(N, W1, A)
    print(f"For the core, the committee W_1 that minimizes satisfaction for N consists of:")
    print(f"- All 8 candidates from C_common: {W1_C_common}")
    print(f"- All 8 candidates from Q_9 and Q_10.")
    print(f"- 4 candidates from the P_i sets (e.g., P_1 and P_2).")
    print(f"The satisfaction s(N,W_1) is {s_N_W1}.\n")


    # 3. Determine s(N, W2) for the EJR committee
    # EJR requires for any group S and integer l, if |S| >= l*n/k, it's not
    # the case that all voters in S approve l common candidates and all have
    # |A(i) intersect W| < l.
    # For S=N={1..8}, l=8: |S|=8 >= 8*10/20=4. They have 8 common candidates (C_common).
    # So, at least one voter in N must have |A(i) intersect W| >= 8.
    # To minimize s(N,W), we want to make the satisfaction values as low as
    # possible while meeting the constraints. This involves an optimization:
    # min s(N,W) = 8*c_w + sum(p_i_w) subject to EJR constraints.
    # This leads to choosing c_w=6, sum(p_i_w)=6, q_9_w=4, q_10_w=4.
    # c_w is the number of candidates from C_common in W.
    
    # We construct W2 to meet these parameters:
    W2_C_common_part = set(range(1, 7))  # 6 candidates from C_common (c_w=6)
    W2_Q = Q_9.union(Q_10) # 8 candidates (q_9_w=4, q_10_w=4)
    # We need sum(p_i_w)=6, with one p_j_w=2 to satisfy the EJR constraint for N.
    # We can pick P[1], P[2], P[3].
    W2_P_part = P[1].union(P[2]).union(P[3]) # 6 candidates
    W2 = W2_C_common_part.union(W2_Q).union(W2_P_part)

    s_N_W2 = calculate_satisfaction(N, W2, A)
    print(f"For EJR, the committee W_2 that minimizes satisfaction for N consists of:")
    print(f"- 6 candidates from C_common.")
    print(f"- All 8 candidates from Q_9 and Q_10.")
    print(f"- 6 candidates from the P_i sets (e.g., P_1, P_2, and P_3).")
    print(f"The satisfaction s(N,W_2) is {s_N_W2}.\n")

    # 4. Calculate and print the ratio
    ratio = Fraction(s_N_W1, s_N_W2)
    print("The ratio is s(N,W_1) / s(N,W_2):")
    print(f"{s_N_W1} / {s_N_W2} = {ratio.numerator}/{ratio.denominator}")

solve_committee_problem()