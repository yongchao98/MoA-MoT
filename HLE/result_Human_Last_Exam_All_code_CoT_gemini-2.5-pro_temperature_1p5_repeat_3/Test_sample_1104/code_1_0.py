def solve_proportionality_puzzle():
    """
    Solves for s1 and s2 based on the rules of PJR and EJR.
    """
    k = 100  # Committee size

    # Part 1: Finding s1 (PJR)
    print("--- Finding s1 for Proportional Justified Representation (PJR) ---")
    
    # To leave voter 1 unsatisfied while satisfying PJR, we must ensure that the group
    # of voters supporting voter 1's ballot does not meet the threshold to demand representation.
    # The set of candidates in voter 1's ballot is C = A(1) = {a, b, c, x}.
    # The number of voters supporting C, |N(C)|, is at least 6.
    min_supporters_of_A1 = 6
    
    # The PJR condition is avoided if |N(A(1))| < s1 / k.
    # This gives the inequality: 6 < s1 / 100.
    # s1 > 600. The smallest integer value is 601.
    s1 = 601
    
    print("The condition to satisfy is: min_supporters_of_A1 < s1 / k")
    print(f"The values are: {min_supporters_of_A1} < {s1} / {k}")
    print(f"This evaluates to: {min_supporters_of_A1 < s1 / k}")
    print(f"The smallest integer s1 is {s1}.\n")
    
    # Part 2: Finding s2 (EJR)
    print("--- Finding s2 for Extended Justified Representation (EJR) ---")
    
    # To leave voter 1 unsatisfied while satisfying EJR, we must ensure that any cohesive
    # group containing only voter 1 does not meet the threshold to demand representation.
    # Voter 1 is in a cohesive group of size 1, V_x = {1}, who all approve candidate 'x'.
    min_size_of_v1_cohesive_group = 1
    
    # The EJR condition for this group is avoided if |V_x| < s2 / k.
    # This gives the inequality: 1 < s2 / 100.
    # s2 > 100. The smallest integer value is 101.
    s2 = 101
    
    print("The condition to satisfy is: min_size_of_v1_cohesive_group < s2 / k")
    print(f"The values are: {min_size_of_v1_cohesive_group} < {s2} / {k}")
    print(f"This evaluates to: {min_size_of_v1_cohesive_group < s2 / k}")
    print(f"The smallest integer s2 is {s2}.\n")
    
    # Final answer
    solution = (s1, s2)
    print(f"The final solution is the pair (s1, s2) = {solution}")
    
    # Required output format
    print(f'<<<{solution}>>>')

solve_proportionality_puzzle()