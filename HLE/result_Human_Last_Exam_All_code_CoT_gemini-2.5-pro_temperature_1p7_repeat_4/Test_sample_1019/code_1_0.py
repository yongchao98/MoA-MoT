import math

def print_hamiltonicity_threshold_formula():
    """
    This function prints the symbolic formula for the d-threshold for Hamiltonicity.
    The problem defines a graph H_n with minimum degree d >= n/2 - eta.
    The d-threshold is the probability p such that H_n U G(n, p) is Hamiltonian
    asymptotically almost surely.

    The final formula is derived based on the worst-case graph H_n, which consists
    of two disjoint cliques of sizes n1 and n2, chosen to meet the minimum degree condition.
    n1 = n/2 - eta + 1
    n2 = n/2 + eta - 1

    The threshold for Hamiltonicity of K_n1 U K_n2 U G(n, p) is determined by the threshold
    for the random bipartite graph between the two cliques to achieve minimum degree 1.
    This threshold is p ~ log(n2)/n1.
    """
    
    # Define the components of the formula as strings
    numerator = "log(n/2 + \u03B7 - 1)"
    denominator = "n/2 - \u03B7 + 1"
    
    # Print the equation
    # We use 'p(n, eta)' to denote that the threshold p is a function of n and eta.
    # The numbers in the equation are 2 and 1. eta is represented by its Unicode character.
    print(f"The d-threshold for Hamiltonicity is given by the formula:")
    print(f"p(n, \u03B7) = {numerator} / ({denominator})")
    print("\nWhere:")
    print("  n is the number of vertices.")
    print("  \u03B7 (eta) is the parameter in the minimum degree d = n/2 - \u03B7.")
    print("  'log' denotes the natural logarithm.")

print_hamiltonicity_threshold_formula()