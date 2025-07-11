def solve_topology_problem():
    """
    Calculates the number of totally bounded group topologies on the integers
    with no nontrivial convergent sequences, under the assumption of first-countability.
    """

    # Step 1: Consider the condition "no nontrivial convergent sequences" for a
    # first-countable topology. This property implies that the topology must be the discrete one.
    # There is only one discrete topology on the integers.
    num_candidate_topologies = 1

    # Step 2: Check if this candidate topology (the discrete one) is "totally bounded" on Z.
    # The discrete topology on an infinite group like the integers is not totally bounded.
    # We use 0 to represent 'False' for this condition check.
    candidate_is_totally_bounded = 0

    # Step 3: Conclude by combining the findings.
    # The number of topologies satisfying both conditions is the number of candidates
    # that pass the final check.
    final_count = num_candidate_topologies * candidate_is_totally_bounded

    # The resulting equation is 1 * 0 = 0.
    # As requested, we print each number in this final equation.
    # The numbers are: the '1' from num_candidate_topologies,
    # the '0' from candidate_is_totally_bounded, and the '0' from final_count.
    print(num_candidate_topologies)
    print(candidate_is_totally_bounded)
    print(final_count)

solve_topology_problem()