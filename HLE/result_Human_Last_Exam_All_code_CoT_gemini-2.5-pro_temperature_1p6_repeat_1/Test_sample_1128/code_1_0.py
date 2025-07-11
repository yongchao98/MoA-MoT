import math
import sys

def solve():
    """
    Calculates the number of "cool strings" of maximal length for n symbols.
    The maximal length of a cool string with n symbols is 3n.
    A string is cool if between any two identical symbols, there are no other
    two identical symbols.
    The number of such strings is given by the formula: n * (n+1)! / 2.
    """
    try:
        # Read the number of symbols, n, from user input.
        # For demonstration purposes, let's use a fixed n, e.g., n=4.
        # To make this code runnable by the user, you can replace the fixed value
        # with n = int(input("Enter the number of symbols (n): "))
        n_str = input("Enter the number of symbols (n): ")
        n = int(n_str)
        if n < 0:
            print("Number of symbols cannot be negative.", file=sys.stderr)
            return

    except (ValueError, IndexError):
        # Default value if input is not provided or invalid
        print("Invalid or no input provided. Using n=4 as an example.", file=sys.stderr)
        n = 4
    
    # Calculate (n+1)!
    if n + 1 < 0:
      # Factorial is not defined for negative numbers.
      # This case is handled by n < 0 check, but as a safeguard.
      fact_val = "undefined"
      result = "undefined"
      print(f"The number of cool strings for n={n} is: {result}")
    elif n + 1 > 20: # math.factorial might be slow/huge for very large numbers
        print("For n > 20, the result is extremely large and may take time to compute.", file=sys.stderr)
        fact_val = math.factorial(n + 1)
        # In Python, integers have arbitrary precision, so division is fine.
        result = (n * fact_val) // 2
    else:
        fact_val = math.factorial(n + 1)
        result = (n * fact_val) // 2

    # Print the equation with the calculated values
    print(f"For n = {n}:")
    print("The formula is: n * (n+1)! / 2")
    print("Step 1: Calculate (n+1)!")
    print(f"({n}+1)! = {n+1}! = {fact_val}")
    print("\nStep 2: Plug values into the formula")
    print(f"Expression: {n} * {fact_val} // 2")
    print(f"Result: {result}")


solve()
<<<36>>>