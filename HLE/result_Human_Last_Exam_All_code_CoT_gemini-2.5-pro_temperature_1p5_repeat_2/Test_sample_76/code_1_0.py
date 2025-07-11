import math

def find_integer_solutions():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    count = 0
    solutions = []
    # Based on analysis, we expect solutions in a relatively small range.
    # A search up to 2000 is more than sufficient.
    for x in range(-8, 2000):
        p_x = x**3 - 16*x**2 - 72*x + 1056
        if p_x >= 0:
            # Check if p_x is a perfect square
            y = math.isqrt(p_x)
            if y * y == p_x:
                count += 1
                solutions.append((x, y))
                
    if solutions:
        print("Found solutions for x:")
        for x, y in solutions:
            # (-4)^3 - 16*(-4)^2 - 72*(-4) + 1056 = -64 - 256 + 288 + 1056 = 1024 = 32^2
            # 4^3 - 16*4^2 - 72*4 + 1056 = 64 - 256 - 288 + 1056 = 576 = 24^2
            # 17^3 - 16*17^2 - 72*17 + 1056 = 4913 - 4624 - 1224 + 1056 = 121 = 11^2
            # Let's print out the full equation.
            x_term = f"({x})^3" if x < 0 else f"{x}^3"
            x2_term = f"16*({x})^2" if x < 0 else f"16*{x}^2"
            x1_term = f"72*({x})" if x < 0 else f"72*{x}"
            val_x3 = x**3
            val_16x2 = 16*x**2
            val_72x = 72*x
            sign_72x = "-" if val_72x >= 0 else "+"
            val_72x = abs(val_72x)
            
            print(f"For x = {x}:")
            print(f"{x_term} - {x2_term} - {x1_term} + 1056 = {val_x3} - {val_16x2} - {72*x} + 1056 = {p_x} = {y}^2")
        print(f"\nThe number of integers x for which the quantity is a perfect square is {count}.")
    else:
        print("No integer solutions found in the checked range.")

find_integer_solutions()