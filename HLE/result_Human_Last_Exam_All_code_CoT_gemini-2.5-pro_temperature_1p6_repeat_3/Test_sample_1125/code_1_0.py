def solve_hat_puzzle():
    """
    This function calculates and explains the largest number of people (N) 
    guaranteed to determine their hat number in the given puzzle.
    """
    total_members = 12
    
    # The core of the strategy is to divide the 12 members into the maximum
    # number of independent, effective groups. An effective group is one that can
    # guarantee at least one member's success.
    
    # Through analysis, the best strategy for a group is to form a "k-cycle with a chord".
    # - A simple 3-person cycle (triangle) is not enough, as the leader can create an
    #   ambiguous situation where no one is certain of their number.
    # - A 4-person cycle is the smallest group where a "chord" (an edge between non-adjacent
    #   members) can be added. This structure is robust against the leader's ambiguity.
    min_group_size_for_guarantee = 4
    
    # By creating the maximum number of these small, robust groups, the team maximizes
    # the total number of guaranteed successes.
    num_groups = total_members // min_group_size_for_guarantee
    
    # Now, consider the worst-case for one such group of 4.
    # The leader, wanting to minimize successes, can make revelations on the cycle
    # that are not all distinct (e.g., the number revealed for pair (A,B) and (B,C)
    # is the same). This is only possible if that number belongs to member B.
    # This guarantees member B knows their number. The leader can then make the
    # other revelations in a way that the rest of the group remains uncertain.
    # So, the guaranteed number of successes per group (the worst-case outcome) is 1.
    guaranteed_success_per_group = 1
    
    # The total number of guaranteed successes (N) is the sum of successes from
    # each independent group.
    max_n = num_groups * guaranteed_success_per_group

    print("The optimal strategy is to divide the 12 members into several independent groups.")
    print(f"The smallest group size that can guarantee a success is {min_group_size_for_guarantee}.")
    print(f"From {total_members} members, a maximum of {num_groups} such groups can be formed.")
    print(f"In the worst-case scenario, each group yields {guaranteed_success_per_group} guaranteed success.")
    print("Therefore, the maximum number of people, N, guaranteed to know their number is the total from all groups.")
    
    print("\nFinal Equation:")
    # The final equation demonstrates the calculation based on the strategy.
    print(f"{num_groups} * {guaranteed_success_per_group} = {max_n}")

solve_hat_puzzle()