def solve_topology_problem():
    """
    This function represents the steps to derive the fundamental group
    and prints the final result.
    """
    # Step 1: Identify the fundamental group of the base shapes.
    # The fundamental group of a torus is Z x Z.
    # The fundamental group of a circle is Z.
    group_Z = "Z"
    direct_product = "x"
    free_product = "*"
    
    # Step 2: Form the group for the torus component.
    torus_group = f"({group_Z} {direct_product} {group_Z})"
    
    # Step 3: The final space is a torus wedged with a circle.
    # Its fundamental group is the free product of the individual groups.
    final_group = f"{torus_group} {free_product} {group_Z}"
    
    # Step 4: Print the final equation as the result.
    # The prompt asks to "output each number in the final equation", which
    # is interpreted here as printing the complete symbolic name of the group.
    print(final_group)

solve_topology_problem()