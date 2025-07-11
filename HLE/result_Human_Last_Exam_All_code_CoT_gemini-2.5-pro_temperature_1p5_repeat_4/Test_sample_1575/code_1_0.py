def solve_reversal_puzzle():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """

    # Step 1 & 2: Define problem parameters based on the analysis.
    # The 100 positions are divided into 5 groups based on index mod 5.
    # Free swaps can arrange elements arbitrarily within each group.
    # Costly adjacent swaps move elements between groups.
    num_elements = 100
    num_groups = 5
    group_size = num_elements // num_groups

    # Step 3 & 4: Determine the required migration of elements between groups.
    # An element k (1-indexed) starts in group (k-1)%5 and must end in group (100-k)%5.
    # This leads to a specific permutation of element sets among the groups.
    # migrations maps start_group -> end_group for the 20 elements in start_group.
    migrations = {
        0: 4,  # Group 0's elements must move to Group 4
        1: 3,  # Group 1's elements must move to Group 3
        2: 2,  # Group 2's elements stay in Group 2
        3: 1,  # Group 3's elements must move to Group 1
        4: 0,  # Group 4's elements must move to Group 0
    }

    # Step 5: Calculate the cost for each migration.
    # The cost is the minimum number of swaps needed.
    # Swaps are only between adjacent groups in the cycle 0-1-2-3-4-0.
    
    # Traffic flow forward (i -> i+1) and backward (i+1 -> i) at each boundary.
    # Boundary 0 is between Group 0 and 1, B1 between G1/G2, etc. B4 between G4/G0.
    flow_fwd = [0] * num_groups
    flow_bwd = [0] * num_groups

    # G1 <-> G3 migration: path G1->G2->G3 (2 hops) is shorter than G1->G0->G4->G3 (3 hops)
    # Flow from G1 to G3 contributes to forward flow at boundaries B1 and B2.
    flow_fwd[1] += group_size # G1 -> G2
    flow_fwd[2] += group_size # G2 -> G3
    # Flow from G3 to G1 contributes to backward flow at boundaries B2 and B1.
    flow_bwd[2] += group_size # G3 -> G2
    flow_bwd[1] += group_size # G2 -> G1

    # G0 <-> G4 migration: path G0->G4 (1 hop) is shorter than G0->G1->...->G4 (4 hops)
    # Flow from G4 to G0 contributes to forward flow at boundary B4.
    flow_fwd[4] += group_size # G4 -> G0
    # Flow from G0 to G4 contributes to backward flow at boundary B4.
    flow_bwd[4] += group_size # G0 -> G4

    # Calculate cost at each boundary.
    # Cost = max(flow in one direction, flow in the other)
    costs = [max(f, b) for f, b in zip(flow_fwd, flow_bwd)]
    
    # Step 6: Sum the costs for the total minimum number of moves.
    total_cost = sum(costs)
    
    # Output the final calculation clearly.
    # The boundaries are B0(G0,G1), B1(G1,G2), B2(G2,G3), B3(G3,G4), B4(G4,G0)
    print("The minimum number of moves is the sum of swaps at each group boundary:")
    print(f"Cost(B0) + Cost(B1) + Cost(B2) + Cost(B3) + Cost(B4) = Total Moves")
    # Using ' + '.join to construct the equation string
    equation_parts = [str(c) for c in costs]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_cost}")

solve_reversal_puzzle()
<<<60>>>