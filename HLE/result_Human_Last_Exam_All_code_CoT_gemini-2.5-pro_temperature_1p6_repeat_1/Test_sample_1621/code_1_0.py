def solve_problem():
    """
    This script determines for how many natural numbers n the given condition holds.
    The solution connects the linear algebra problem to a theorem in algebraic topology.
    """
    
    print("The problem is to find the number of natural numbers 'n' for which there exist n real n-by-n matrices")
    print("A_1, ..., A_n such that for any non-zero vector x, the vectors A_1*x, ..., A_n*x are linearly independent.\n")
    
    print("Step 1: The condition on the matrices can be shown to imply that the (n-1)-dimensional sphere, S^(n-1), must be a 'parallelizable manifold'.")
    print("This means it is possible to choose n-1 tangent vector fields at each point of the sphere that are linearly independent.\n")

    print("Step 2: A famous theorem in algebraic topology, proved by J.F. Adams in 1962, gives a complete list of all parallelizable spheres.")
    print("It states that a sphere S^k is parallelizable only if its dimension k is 0, 1, 3, or 7.\n")
    
    parallelizable_sphere_dimensions = [0, 1, 3, 7]
    possible_values_of_n = []

    print("Step 3: In our problem, the dimension of the sphere is k = n-1.")
    print("We can find the possible values for 'n' by solving the equation n - 1 = k for each possible k:")

    for k in parallelizable_sphere_dimensions:
        n = k + 1
        possible_values_of_n.append(n)
        # This part outputs each number in the final set of equations
        print(f"  For k = {k}:  n - 1 = {k}  =>  n = {n}")

    print("\nStep 4: It can be proven that these values of n are indeed possible.")
    print("The required matrices can be constructed using real division algebras, which exist only in these dimensions (1, 2, 4, 8).")
    
    print("\nConclusion:")
    print(f"The only natural numbers n for which such matrices exist are {possible_values_of_n}.")
    print(f"There are {len(possible_values_of_n)} such numbers.")

solve_problem()