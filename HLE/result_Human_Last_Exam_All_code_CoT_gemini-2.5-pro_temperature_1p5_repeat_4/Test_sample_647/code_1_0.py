import math

def solve():
    """
    Solves the committee election problem by analyzing the properties of
    committees in the core and those satisfying EJR to find the ones
    that minimize the satisfaction for group N, and then computes their ratio.
    """

    # --- Problem Definition ---
    # Number of voters
    n = 10
    # Committee size
    k = 20
    # Group of interest is N = {1, ..., 8}

    # Let W be a committee. We define:
    # x = |W intersect {1,...,8}|
    # z = |W intersect {9,...,24}|
    # w = |W intersect {25,...,32}|
    # The satisfaction for group N is s(N,W) = 8*x + z.
    # The committee size is x + z + w = 20.

    # --- Part 1: Finding s(N, W1) for the Core ---

    # A committee W is in the core if it is stable against deviations by groups of voters.
    # It can be shown that any committee W that does not contain all common candidates {1,...,8}
    # can be blocked by a large subgroup of N={1,...,8}. All these voters would prefer
    # to swap a less-liked candidate in W for one of the missing common candidates.
    # Therefore, for any committee W to be in the core, it must contain all 8 common candidates.
    x1 = 8

    # For W1, a core committee with the lowest satisfaction for N:
    # s(N, W1) = 8*x1 + z1 = 8*8 + z1 = 64 + z1.
    # To minimize this value, we must choose the rest of the committee to minimize z1.
    # The committee size is fixed: x1 + z1 + w1 = 20.
    # With x1 = 8, we have 8 + z1 + w1 = 20, which simplifies to z1 + w1 = 12.
    # To minimize z1, we must maximize w1 = |W1 intersect {25,...,32}|.
    # The set {25,...,32} contains 8 candidates, so the maximum possible value for w1 is 8.
    
    w1_max = 8
    w1 = w1_max
    z1 = 12 - w1
    s_N_W1 = 8 * x1 + z1

    print("--- Analysis for W1 (Core) ---")
    print("For any committee W in the core, it must contain all common candidates {1,...,8}.")
    print(f"This means x = |W intersect {{1,...,8}}| = {x1}.")
    print(f"The satisfaction for group N is s(N,W) = 8*x + z = 8*{x1} + z = {8*x1} + z.")
    print("To minimize satisfaction, we must minimize z = |W intersect {9,...,24}|.")
    print("Given x + z + w = 20 and x = 8, we have z + w = 12.")
    print(f"Minimizing z requires maximizing w = |W intersect {{25,...,32}}|.")
    print(f"The maximum possible w is {w1_max}, as there are {w1_max} candidates in {{25,...,32}}.")
    print(f"So, for W1, we have x1 = {x1}, w1 = {w1}, which implies z1 = 12 - {w1} = {z1}.")
    print(f"The lowest satisfaction for N with a core committee W1 is s(N,W1) = 8 * {x1} + {z1} = {s_N_W1}.")
    print("-" * 20)

    # --- Part 2: Finding s(N, W2) for EJR ---

    # A committee W satisfies Extended Justified Representation (EJR) if any sufficiently large
    # and cohesive group of voters gets a proportional number of approved candidates in the committee.
    # Applying the EJR definition to this problem's profile yields several constraints on W.
    # Notably, for each voter i from 1 to 10, |W intersect A(i)| must be at least 2.

    # We want to find an EJR committee W2 that minimizes s(N, W2) = 8*x2 + z2.
    # To minimize this, we should try to make x2 as small as possible. Let's test if x2=0 is possible.
    
    x2_trial = 0
    # If x2 = 0, the EJR constraint for voters i=1..8, which is |W intersect A(i)| >= 2,
    # becomes |W intersect {c_i1, c_i2}| >= 2. Since each pair {c_i1, c_i2} has exactly 2 candidates,
    # W must contain both of them for all i=1..8. This means W must contain all candidates from {9,...,24}.
    z2_forced_by_x0 = 16
    
    # If x2=0 and z2=16, then from x2+z2+w2=20, we get w2 = 20 - 0 - 16 = 4.
    w2_forced_by_x0 = 4
    
    # A committee with w2=4 can satisfy the other EJR constraints (e.g., |W intersect {25..28}| >= 2).
    # Thus, an EJR committee with x2=0 is possible. This choice minimizes 8*x+z.
    
    x2 = 0
    z2 = 16
    s_N_W2 = 8 * x2 + z2

    print("--- Analysis for W2 (EJR) ---")
    print("EJR implies several constraints, including |W intersect A(i)| >= 2 for each voter i=1..10.")
    print(f"To find the EJR committee W2 with the lowest satisfaction s(N,W) = 8*x + z, we should minimize x.")
    print(f"We test if a committee with x = {x2_trial} can satisfy all EJR constraints.")
    print(f"If x = {x2_trial}, EJR constraints for voters i=1..8 force W to include all their unique candidates.")
    print(f"This means W must contain all of {{9,...,24}}, so z = {z2_forced_by_x0}.")
    print(f"From x+z+w=20, this implies w = {w2_forced_by_x0}, which is compatible with all constraints.")
    print(f"So, the EJR committee W2 with the minimum satisfaction for N has x2 = {x2} and z2 = {z2}.")
    print(f"The lowest satisfaction for N with an EJR committee W2 is s(N,W2) = 8 * {x2} + {z2} = {s_N_W2}.")
    print("-" * 20)

    # --- Final Calculation ---
    ratio = s_N_W1 / s_N_W2
    print("--- Final Result ---")
    print(f"The satisfaction for N with W1 is s(N,W1) = {s_N_W1}.")
    print(f"The satisfaction for N with W2 is s(N,W2) = {s_N_W2}.")
    print(f"The final equation is s(N,W1) / s(N,W2) = {s_N_W1} / {s_N_W2}.")
    print(f"The ratio is {ratio}")
    
    return ratio

final_ratio = solve()
print(f"<<<{final_ratio}>>>")