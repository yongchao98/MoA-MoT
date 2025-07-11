import math

def solve_cake_cutting_bound():
    """
    Calculates the n-dependent part of the query complexity bound for
    the 4-agent connected ε-envy-free cake-cutting problem.
    """
    # The number of agents in the problem.
    n = 4

    # The query complexity for finding a connected ε-envy-free allocation for n
    # agents is O(n*log(n)/ε). We calculate the n-dependent factor.
    # In complexity analysis, log is conventionally log base 2.
    log_n = math.log2(n)

    # The resulting upper bound constant O for n=4.
    result = n * log_n

    print("For the envy-free cake-cutting problem with n agents, the tightest")
    print("known upper bound on query complexity for a connected ε-envy-free")
    print("allocation is given by the formula O(n * log(n) / ε).")
    print("\nWe calculate the n-dependent factor for n = 4, using log base 2:")
    
    # Print the equation with all the numbers involved.
    print(f"\nFinal Equation: {int(n)} * log2({int(n)}) = {int(n)} * {int(log_n)} = {int(result)}")

solve_cake_cutting_bound()
