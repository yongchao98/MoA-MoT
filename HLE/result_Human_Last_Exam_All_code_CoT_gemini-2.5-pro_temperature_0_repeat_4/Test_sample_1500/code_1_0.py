def solve_lie_group_questions():
    """
    Provides answers and demonstrations for the questions on coadjoint orbits.
    """
    
    # Print the answers to the three questions first.
    print("Answers:")
    print("(a) True")
    print("(b) No")
    print("(c) No")
    print("\n" + "="*40 + "\n")
    
    # --- Demonstration for part (b) ---
    print("Demonstration for (b): Is b_2(O_lambda) for SU(n) always n-1?")
    print("The formula for the second Betti number is b_2 = k - 1, where k is the number of distinct eigenvalues of lambda.")
    
    n = 4
    # Consider a singular element lambda for SU(4) where eigenvalues have multiplicities [2, 2].
    # For example, lambda could be diag(i, i, -i, -i).
    # The number of distinct eigenvalues is k=2.
    partition = [2, 2]
    k = len(partition)
    
    # Calculate b_2 using the formula
    b2_calculated = k - 1
    
    # The value claimed in the question
    n_minus_1 = n - 1
    
    print(f"\nLet's take G = SU({n}).")
    print(f"Consider a coadjoint orbit O_lambda where lambda has eigenvalue multiplicities given by the partition {partition}.")
    print(f"This means there are k = {k} distinct eigenvalues.")
    print("The second Betti number is given by the equation: b_2 = k - 1.")
    print(f"Substituting k = {k}, we get: b_2 = {k} - 1 = {b2_calculated}.")
    print(f"The question suggests the value should be n - 1 = {n} - 1 = {n_minus_1}.")
    print(f"Since {b2_calculated} != {n_minus_1}, the statement is false.\n")
    
    print("="*40 + "\n")

    # --- Demonstration for part (c) ---
    print("Demonstration for (c): Must H_G*(O_lambda) be isomorphic to a GKM graph cohomology?")
    print("This requires the manifold O_lambda to satisfy the GKM conditions.")
    print("A key GKM condition is that the isotropy weights at each T-fixed point must be pairwise linearly independent.")
    
    print("\nLet's take G = SU(3) and the full flag manifold O_lambda = SU(3)/T^2 (where lambda is regular).")
    print("At the identity fixed point, the isotropy weights are the roots of the Lie algebra su(3).")
    
    # The roots of su(3) can be represented as vectors e_i - e_j.
    # Let's represent them symbolically.
    roots = ["(e1-e2)", "(e2-e1)", "(e1-e3)", "(e3-e1)", "(e2-e3)", "(e3-e2)"]
    
    print("The set of roots for su(3) is: {", ", ".join(roots), "}")
    
    root1 = "(e1-e2)"
    root2 = "(e2-e1)"
    
    print(f"This set contains, for example, the root {root1} and the root {root2}.")
    print(f"Note that {root2} = -{root1}.")
    print("Since one root is a scalar multiple of another, they are linearly dependent.")
    print("The GKM condition is therefore not satisfied.")
    print("Thus, the statement that the equivariant cohomology *must* be isomorphic to that of a GKM graph is false.")

# Execute the function to print the solution.
solve_lie_group_questions()