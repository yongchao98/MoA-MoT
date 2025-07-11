import math

def solve_envy_free_bound():
    """
    Calculates the upper bound O for a connected epsilon-envy-free
    allocation for four agents based on a known complexity formula.

    The query complexity for n agents can be bounded by O(n^n * poly(1/epsilon)).
    To find a concrete numerical bound 'O', we evaluate the n^n component for n=4.
    """
    
    # Number of agents
    n = 4
    
    # Calculate the upper bound based on the n^n formula
    upper_bound = int(math.pow(n, n))
    
    # The final equation requires printing each number
    base = n
    exponent = n
    
    # Print the equation and the result
    print(f"The upper bound O can be estimated using the formula n^n.")
    print(f"For n = 4 agents:")
    print(f"{base} ^ {exponent} = {upper_bound}")

solve_envy_free_bound()