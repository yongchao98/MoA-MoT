from sympy import Function, Symbol, Integer

def print_obstruction_groups():
    """
    This function symbolically represents and prints the homotopy-theoretic
    obstruction groups for the given problem.
    """
    # Define symbolic representations for mathematical objects
    H_tilde = Function('H_tilde') # Reduced homology group
    pi = Function('pi')         # Homotopy group
    SO = Function('SO')         # Special Orthogonal group
    Hom = Function('Hom')       # Group of homomorphisms
    
    # Define symbols for variables in the problem
    n = Symbol('n')
    k = Symbol('k')
    X = Symbol('X')

    # The first obstruction group is the fundamental group of SO(2k)
    group1 = pi(1, SO(2 * k))
    
    # The second obstruction group is related to the homology of X and
    # a higher homotopy group of SO(2k).
    # H_tilde_{n-1}(X) is the (n-1)-th reduced homology group of X.
    # Since X is a homology (n-1)-sphere, this group is isomorphic to Z.
    # Hom(Z, G) is isomorphic to G. We write it in the more descriptive
    # form involving the homology of X as requested.
    homology_group_X = H_tilde(n - 1, X)
    homotopy_group_SO = pi(n + 1, SO(2 * k))
    group2 = Hom(homology_group_X, homotopy_group_SO)

    # Print the results
    print("The homotopy-theoretic obstructions for the two paths to be homotopic are elements")
    print("of a group which is an extension of the following two groups:")
    print("\nObstruction Group 1:")
    print(f"  {group1}")
    print("\nObstruction Group 2:")
    print(f"  {group2}")
    
    print("\nNote: Since X is a homology (n-1)-sphere, "
          f"{homology_group_X} is isomorphic to the integers Z.")
    print("Therefore, Obstruction Group 2 is isomorphic to "
          f"{homotopy_group_SO}.")
    
# Execute the function to display the answer
print_obstruction_groups()
