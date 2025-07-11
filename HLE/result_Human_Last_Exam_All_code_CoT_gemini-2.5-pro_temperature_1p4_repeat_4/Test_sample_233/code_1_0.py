import sys

def solve():
    """
    Calculates the smallest positive integer g for the given topology problem.
    """
    # The given surface is Sigma.
    # Genus of Sigma.
    g_sigma = 10
    # Number of boundary components of Sigma.
    b_sigma = 1

    # The problem asks for the smallest genus g of a closed surface Sigma'
    # that contains Sigma, regardless of how Sigma is embedded.

    # A general construction for such a containing surface is the "double" of Sigma.
    # The double of a surface S, denoted D(S), is a closed surface whose genus
    # is given by the formula: g(D(S)) = 2 * g(S) + b(S) - 1,
    # where g(S) is the genus of S and b(S) is the number of its boundary components.

    # This construction provides an upper bound for g. It can be shown in topology
    # that this bound is "sharp", meaning there exist "worst-case" embeddings of Sigma
    # that require a containing surface of this exact genus.
    # Therefore, the smallest g that works for all cases is the genus of the double.

    g_final = 2 * g_sigma + b_sigma - 1
    
    # Print the explanation and the final equation.
    print("The problem requires finding the smallest genus g for a closed surface that can contain any given surface Σ (genus 10, 1 boundary).")
    print("This can be determined by calculating the genus of the topological 'double' of Σ.")
    print("The formula for the genus of a double (g') from a surface with genus g and b boundaries is: g' = 2*g + b - 1.")
    print("\nFor our surface Σ:")
    print(f"g = {g_sigma}")
    print(f"b = {b_sigma}")
    print("\nPlugging these values into the formula:")
    print(f"g' = 2 * {g_sigma} + {b_sigma} - 1 = {g_final}")
    print(f"\nThus, the smallest positive integer g is {g_final}.")

solve()