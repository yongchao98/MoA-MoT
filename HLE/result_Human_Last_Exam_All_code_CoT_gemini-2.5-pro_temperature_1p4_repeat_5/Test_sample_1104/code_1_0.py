def solve_election_proportionality():
    """
    This script calculates s1 and s2 based on the definitions of PJR and EJR.
    It explains the logic and shows the numbers used in the final equations.
    """

    # --- Problem Parameters ---
    # Committee size
    k = 100

    # --- Part 1: Proportional Justified Representation (PJR) ---
    print("### Solving for s_1 (PJR) ###\n")
    print("A committee satisfies PJR if for any cohesive group of unsatisfied voters N', its size l must be less than n/k.")
    print("To find the smallest n (s1) that allows voter 1 to be unsatisfied, we must construct a committee W")
    print("such that the size of the smallest cohesive group of unsatisfied voters containing voter 1, |U_min|, is minimized.\n")

    print("We are given that voters 1-6 all approve {a,b,c}, forming a cohesive group.")
    print("To make voter 1 unsatisfied, W cannot contain {a,b,c,x}.")
    print("To minimize the number of unsatisfied voters in this group, W should satisfy voters 2-6.")
    print("This is achieved by including 'y' (satisfies voters 2,3) and 'z' (satisfies voters 4,5,6) in W.")
    min_U_size = 1
    print(f"This leaves only voter 1 unsatisfied, so the minimum group size is |U_min| = {min_U_size}.\n")

    print("The PJR condition that must be satisfied is: |U_min| < n / k")
    print(f"Substituting the values, we get the inequality: {min_U_size} < n / {k}")
    print(f"This simplifies to n > {min_U_size} * {k}, which is n > {min_U_size * k}.")
    s1 = min_U_size * k + 1
    print("The final equation for s1 is derived from finding the smallest integer n satisfying this.")
    print(f"s1 = {min_U_size} * {k} + 1 = {s1}\n\n")

    # --- Part 2: Extended Justified Representation (EJR) ---
    print("### Solving for s_2 (EJR) ###\n")
    print("A committee violates EJR if there exists a j-cohesive group N' of size l >= n/k where all voters i in N' have fewer than j of their approved candidates in W.\n")
    
    print("Consider the group N' of the first 6 voters. Its size is l = 6.")
    l = 6
    print("This group is 3-cohesive because they all approve of the 3 candidates {a,b,c}. So, j = 3.")
    j = 3

    print(f"\nAn EJR violation with j={j} occurs if two conditions are met:")
    print(f"1. Size condition: l >= n/k  =>  {l} >= n/{k}  =>  n <= {l * k}")
    print(f"2. Satisfaction condition: For all i in N', |A(i) intersect W| < {j}.\n")
    
    print("We check the satisfaction condition, assuming voter 1 is unsatisfied:")
    print("  - For voter 1: |A(1) intersect W| = 0, which is < 3.")
    print("  - For voters 2,3: At most, 'y' can be in W, so |A(i) intersect W| <= 1, which is < 3.")
    print("  - For voters 4,5,6: At most, 'z' can be in W, so |A(i) intersect W| <= 1, which is < 3.")
    print("This means the satisfaction condition always holds for this group with j=3.\n")
    
    print("Thus, an EJR violation is triggered if and only if the size condition n <= 600 holds.")
    print("To satisfy EJR, we must ensure n > 600.")
    print(f"The inequality to avoid the violation is: l < n / k")
    print(f"Substituting the values: {l} < n / {k}")
    print(f"This simplifies to n > {l} * {k}, which is n > {l * k}.")
    s2 = l * k + 1
    print("The final equation for s2 is derived from finding the smallest integer n satisfying this.")
    print(f"s2 = {l} * {k} + 1 = {s2}\n")
    
    final_answer = (s1, s2)
    print(f"The solution is the pair (s1, s2) = {final_answer}.")

solve_election_proportionality()