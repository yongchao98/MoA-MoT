import sys

# This problem is a classic example from algebraic topology.
# The code will not perform a complex calculation but will print the logical steps to derive the answer.

def solve_topology_problem():
    """
    Derives the fundamental group of the described topological space.
    """
    
    # Step 1: Define the fundamental group of a single pair of pants.
    # A pair of pants (a 3-holed sphere) deformation retracts to a wedge of 2 circles.
    # Its fundamental group is the free group on 2 generators.
    group_pants = "Z * Z"
    num_generators_per_pant = 2
    
    # Let the generators for the first pair of pants be a, b.
    # Let the generators for the second pair of pants be c, d.
    
    # Step 2: Combine the two pairs of pants at the waistbands.
    # Identifying the waistbands to a single point creates a wedge sum of the two spaces.
    # The fundamental group of a wedge sum is the free product of the individual groups.
    num_initial_generators = num_generators_per_pant * 2
    initial_group = " * ".join(["Z"] * num_initial_generators)
    
    print(f"Step 1: The fundamental group of a single pair of pants can be represented by {num_generators_per_pant} generators. For two pants, we initially have a total of {num_initial_generators} generators.")
    print(f"The group is the free product of the two, resulting in a group with {num_initial_generators} generators: {initial_group}.")
    print("Let the generators be a, b, c, d.\n")
    
    # Step 3: Apply the "sewing" relations.
    # Sewing the corresponding leg openings together introduces relations between the generators.
    # Sewing leg 'a' to 'c' implies the relation a = c.
    # Sewing leg 'b' to 'd' implies the relation b = d.
    num_relations = 2
    
    print("Step 2: Sewing the legs together introduces relations.")
    print("Relation 1: a = c")
    print("Relation 2: b = d\n")

    # Step 4: Calculate the final group.
    # We started with 4 generators and introduced 2 independent relations.
    # This allows us to eliminate 2 generators.
    num_final_generators = num_initial_generators - num_relations
    final_group = " * ".join(["Z"] * num_final_generators)

    print(f"Step 3: The initial group was <a, b, c, d>. With the relations a=c and b=d, we can simplify this to <a, b>.")
    print(f"The final group is the free group on {num_final_generators} generators.")
    sys.stdout.write("The fundamental group is: ")
    
    # Output each Z in the final equation.
    for i in range(num_final_generators):
        sys.stdout.write("Z")
        if i < num_final_generators - 1:
            sys.stdout.write(" * ")
    print()


solve_topology_problem()
