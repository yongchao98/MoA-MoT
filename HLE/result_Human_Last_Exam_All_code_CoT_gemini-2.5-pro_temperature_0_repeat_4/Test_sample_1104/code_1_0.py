import math

def solve_proportionality_problem():
    """
    This script calculates the smallest preference profile sizes (s1, s2)
    for PJR and EJR based on the problem description.
    """
    k = 100  # Committee size

    # --- Part 1: Calculation for s1 (Proportional Justified Representation) ---

    # For PJR, a cohesive group of size l >= n/k forces one of their
    # approved candidates into the committee.
    # To leave voter 1 unsatisfied, their cohesive group must not meet this threshold.
    # The cohesive group for A(1) has size l=1.
    l1 = 1

    # The condition is l < s1 / k, which prevents PJR from being triggered for voter 1.
    # s1 > l * k
    s1 = math.floor(k / l1) + 1

    print("Calculation for s1 (PJR):")
    print("The condition to leave voter 1 unsatisfied is that their cohesive group (size l=1) does not meet the PJR threshold.")
    print(f"l < s1 / k")
    print(f"{l1} < s1 / {k}")
    print(f"s1 > {k * l1}")
    print(f"The smallest integer s1 is {s1}.")
    print("-" * 30)

    # --- Part 2: Calculation for s2 (Extended Justified Representation) ---

    # For EJR, a group of size l with a common approval set C, where |C| >= l*k/n,
    # forces a candidate from C into the committee.
    # To leave voter 1 unsatisfied, any group containing voter 1 must not meet this threshold.
    # This implies n < l*k / |C|. We need to find the tightest constraint.
    
    initial_ballots = [
        {'a', 'b', 'c', 'x'},
        {'a', 'b', 'c', 'y'},
        {'a', 'b', 'c', 'y'},
        {'a', 'b', 'c', 'z'},
        {'a', 'b', 'c', 'z'},
        {'a', 'b', 'c', 'z'},
    ]

    # The tightest constraint comes from the group N'={voter 1}.
    l_ejr = 1
    # The size of the intersection for this group is |A(1)|.
    C_size_ejr = len(initial_ballots[0])

    # The condition is s2 < l * k / |C|.
    min_bound = (l_ejr * k) / C_size_ejr

    # The smallest possible profile size is the initial size.
    s2 = len(initial_ballots)
    # We must check if this smallest size is valid. It is, since s2 < min_bound.

    print("Calculation for s2 (EJR):")
    print("The condition to leave voter 1 unsatisfied is that for any group N' including voter 1, the EJR condition is not met.")
    print("This means s2 < l * k / |C| for all such groups.")
    print("The tightest constraint is from the group N'={voter 1} with l=1 and |C|=4.")
    print(f"s2 < l * k / |C|")
    print(f"s2 < {l_ejr} * {k} / {C_size_ejr}")
    print(f"s2 < {min_bound}")
    print(f"The profile must extend the initial 6 ballots, so the smallest possible size is {s2}.")
    print(f"Since {s2} is less than {min_bound}, a valid construction is possible. Thus, the smallest size s2 is {s2}.")
    print("-" * 30)

    print(f"The final answer is the pair (s1, s2) = ({s1}, {s2}).")
    
    # Final answer in the required format
    print(f"<<<({s1}, {s2})>>>")

solve_proportionality_problem()