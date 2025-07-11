import math

def search_smallest_r():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    (p+1)^2 + (q+1)^2 = (r+1)^2 + 1 for some integers p, q > 1.
    """
    # We are looking for the smallest r, so we iterate on r starting from 2.
    r = 2
    while True:
        # Let z = r+1. We need to solve x^2+y^2 = z^2+1 for x=p+1, y=q+1.
        # So x,y >= 3.
        z = r + 1
        target = z**2 + 1

        # Search for integer pairs (x, y) such that x^2+y^2 = target.
        # To minimize r, we should search for p,q.
        # p starts from 2, so x starts from 3.
        found = False
        for p in range(2, r + 2): # Heuristic limit for p
            x = p + 1
            y_squared = target - x**2
            if y_squared > 0:
                y = math.isqrt(y_squared)
                if y * y == y_squared:
                    # y=q+1, so q=y-1. We need q > 1, so y > 2.
                    if y > 2:
                        q = y - 1
                        print(f"Found a solution for r={r}: p={p}, q={q}.")
                        print(f"Verification: ({p}+1)^2 + ({q}+1)^2 = {x**2} + {y**2} = {x**2 + y**2}")
                        print(f"And ({r}+1)^2 + 1 = {z**2} + 1 = {z**2 + 1}")
                        print(f"The equation holds: {x**2} + {y**2} = {z**2 + 1}")
                        print(f"The smallest such value for r is {r}.")
                        return r
        r += 1

# Execute the search and print the result.
smallest_r = search_smallest_r()
print(f"\nFinal Answer: The smallest possible value of r is {smallest_r}.")