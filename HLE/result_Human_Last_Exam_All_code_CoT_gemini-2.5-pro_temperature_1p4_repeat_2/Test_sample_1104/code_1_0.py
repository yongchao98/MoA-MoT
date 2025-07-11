import math

def solve_proportionality_riddles():
    """
    Solves for the smallest preference profile sizes s1 and s2.
    """
    
    # Committee size
    k = 100

    # From the problem description, we have a group of 6 voters, {1, 2, 3, 4, 5, 6},
    # that all approve of candidate 'a' (and also 'b' and 'c').
    # Let's call the size of this group 'g'.
    g = 6

    # --- Part 1: Proportional Justified Representation (PJR) ---
    # We want voter 1 to be unsatisfied, meaning W must not contain {a,b,c,x}.
    # Let's focus on candidate 'a'. For 'a' not to be in W, the group supporting 'a'
    # must be smaller than the PJR representation threshold n/k.
    # The inequality is g < n/k.
    # 6 < s1 / 100  => s1 > 600
    # The smallest integer s1 is 601.
    s1 = math.floor(g * k) + 1
    
    print("Step 1: Calculating s1 for PJR")
    print(f"The size of the group unanimously supporting candidate 'a' is {g}.")
    print(f"The committee size k is {k}.")
    print("For PJR to allow 'a' to be excluded from the committee, the group size must be less than the threshold n/k.")
    print(f"This gives the inequality: {g} < s1 / {k}")
    print(f"Solving for s1: s1 > {g * k}")
    print(f"The smallest integer value for s1 is therefore {s1}.")
    print("-" * 20)

    # --- Part 2: Extended Justified Representation (EJR) ---
    # For EJR, we test for a group V' of size |V'| >= l*n/k that supports l candidates
    # who are not in W.
    # Let's consider the simplest case, l=1, and candidate 'a'.
    # We want 'a' not in W. The group supporting 'a' has size g=6.
    # For EJR not to be violated, the group size must be smaller than the EJR threshold l*n/k.
    # The inequality is g < l * n / k, with l=1.
    # 6 < 1 * s2 / 100 => s2 > 600
    # The smallest integer s2 is 601.
    l = 1
    s2 = math.floor(g * k / l) + 1

    print("Step 2: Calculating s2 for EJR")
    print(f"The size of the group unanimously supporting candidate 'a' is {g}.")
    print("We test the EJR condition for l=1 and candidate set C={'a'}.")
    print("For EJR to allow 'a' to be excluded, the group size must be less than the threshold l*n/k.")
    print(f"This gives the inequality: {g} < {l} * s2 / {k}")
    print(f"Solving for s2: s2 > {g * k / l}")
    print(f"The smallest integer value for s2 is therefore {s2}.")
    print("-" * 20)
    
    print(f"The final solution is the pair (s1, s2).")
    print(f"s1 = {s1}")
    print(f"s2 = {s2}")

solve_proportionality_riddles()

# The final answer is the pair (s1, s2)
print("<<<({}, {})>>>".format(601, 601))