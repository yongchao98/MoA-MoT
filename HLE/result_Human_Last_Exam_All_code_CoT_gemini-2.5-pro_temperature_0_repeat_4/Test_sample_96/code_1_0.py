import math

def solve_artin_group_problem():
    """
    Calculates the number of torsion elements of order 10 with minimal positive word length
    in the quotient group A/Z for the Artin group A of type E8.
    """
    # Step 1: Identify Key Parameters for the E8 Artin group.
    # The E8 group is defined by 8 standard generators.
    n = 8
    # The Coxeter number for E8 is 30.
    h = 30
    
    print(f"The Artin group of type E8 has n = {n} generators.")
    print(f"The Coxeter number for E8 is h = {h}.")
    print("-" * 20)

    # Step 2: Determine the minimal word length.
    # The desired order of the torsion element is 10.
    d = 10
    print(f"We are looking for torsion elements of order d = {d}.")

    # The order d is related to an integer k by the formula d = h / gcd(h, k).
    # We need to find the smallest positive integer k that satisfies this.
    # 10 = 30 / gcd(30, k)  =>  gcd(30, k) = 3.
    # The smallest positive integer k such that gcd(30, k) = 3 is k itself.
    k = 3
    print(f"The smallest integer k satisfying 10 = 30 / gcd(30, k) is k = {k}.")

    # The minimal length L of a positive word for such an element is L = k * n.
    min_length = k * n
    print(f"The minimal word length is L = k * n = {k} * {n} = {min_length}.")
    print("-" * 20)

    # Step 3 & 4: Identify and count the elements.
    # The elements are of the form c^k, where c is a Coxeter element.
    # We need to count the number of distinct Coxeter elements in the group.
    print("The elements are of the form c^3, where c is a Coxeter element.")
    print("The number of such elements is the number of distinct Coxeter elements.")
    
    # Step 5: The number of distinct Coxeter elements is the number of acyclic orientations
    # of the E8 Coxeter-Dynkin diagram.
    print("This count is equal to the number of acyclic orientations of the E8 Dynkin diagram.")

    # Step 6: The E8 diagram is a tree with n=8 vertices. A tree with n vertices has n-1 edges.
    num_vertices = n
    num_edges = num_vertices - 1
    print(f"The E8 diagram is a tree with {num_vertices} vertices and {num_edges} edges.")

    # For a tree, any orientation of its edges is acyclic.
    # Each of the (n-1) edges can be oriented in 2 ways.
    num_orientations = 2**num_edges
    
    print(f"The number of acyclic orientations is 2^(number of edges).")
    print(f"Final calculation: 2^{num_edges} = {num_orientations}")
    
    print("-" * 20)
    print(f"The number of such torsion elements is {num_orientations}.")

solve_artin_group_problem()