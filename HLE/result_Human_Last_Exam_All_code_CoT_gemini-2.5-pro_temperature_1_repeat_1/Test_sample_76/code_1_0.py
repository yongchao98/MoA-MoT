import math

def solve():
    """
    This function finds all integer values of x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # Define the ranges to check based on analysis
    # Range 1: [-8, 7]
    # Range 2: [17, 62]
    # Range 3: Perfect squares k^2 for k from ceil(sqrt(17)) to 68.
    # The code below combines these checks efficiently.
    
    checked_numbers = set()
    
    # Check the first interval
    for x in range(-8, 8):
        if x in checked_numbers:
            continue
        val = x**3 - 16*x**2 - 72*x + 1056
        if val >= 0:
            s = int(math.sqrt(val))
            if s*s == val:
                solutions.append(x)
        checked_numbers.add(x)
    
    # Check the second interval and the range of squares
    # The upper limit is 68^2 = 4624, since for x >= 69^2, there are no solutions.
    # The analysis showed that for x>=63, we only need to check perfect squares.
    # The loop until 4625 covers all remaining candidates.
    
    for x in range(17, 4625):
        if x in checked_numbers:
            continue
        
        # For x>=63, only check if x is a perfect square
        if x >= 63:
            sqrt_x = int(math.sqrt(x))
            if sqrt_x * sqrt_x != x:
                continue

        val = x**3 - 16*x**2 - 72*x + 1056
        if val >= 0:
            s = int(math.sqrt(val))
            if s*s == val:
                solutions.append(x)
        checked_numbers.add(x)

    print("The integer values of x for which the expression is a perfect square are:")
    for sol in sorted(solutions):
        val = sol**3 - 16*sol**2 - 72*sol + 1056
        s = int(math.sqrt(val))
        print(f"x = {sol}, expression = {val} = {s}^2")
    
    print(f"\nTotal number of integers is: {len(solutions)}")

solve()