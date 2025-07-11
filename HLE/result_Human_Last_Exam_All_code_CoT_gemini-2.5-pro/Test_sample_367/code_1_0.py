import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def gcd_four(a, b, c, d):
    """Computes the greatest common divisor of four numbers."""
    return gcd(gcd(a, b), gcd(c, d))

def solve():
    """
    Finds the number of distinct parallelograms based on the given restrictions.
    """
    solutions = set()
    limit = 100

    # Iterate through possible values for u, v, w
    # v loop
    for v in range(1, int(limit**0.5) + 1):
        # u loop
        for u in range(v + 1, limit):
            # w loop
            for w in range(1, limit):
                # Condition: u, v, w must be pairwise coprime
                if gcd(u, v) != 1 or gcd(u, w) != 1 or gcd(v, w) != 1:
                    continue
                
                # Condition: not a rectangle (v != w)
                if v == w:
                    continue

                # Condition: positive diagonals (uv > w^2)
                if u * v <= w**2:
                    continue

                # Parametrization formulas
                a_prime = u * (v**2 + w**2)
                b_prime = v * (u**2 + w**2)
                
                # Check if a_prime + b_prime is already too large to speed up
                if a_prime + b_prime > 4 * limit : # Heuristic bound
                    if u > v*5 and u > w*5: # if u is much larger, it will only grow
                        break

                d1_prime = abs((u - v) * (u * v + w**2))
                d2_prime = (u + v) * (u * v - w**2)

                # Find primitive parallelogram by dividing by gcd of all components
                common_divisor = gcd_four(a_prime, b_prime, d1_prime, d2_prime)
                
                a = a_prime // common_divisor
                b = b_prime // common_divisor

                # Ensure a < b
                if a > b:
                    a, b = b, a

                # Check final conditions
                # 1. a and b are coprime
                if gcd(a, b) != 1:
                    continue
                
                # 2. 2a < a+b < 100 => a < b and a+b < 100
                if a + b < limit:
                    solutions.add((a, b))

    sorted_solutions = sorted(list(solutions))
    
    print("Found distinct parallelograms for the following (a, b) pairs:")
    for a, b in sorted_solutions:
        print(f"a = {a}, b = {b}, a+b = {a+b}")
    
    print(f"\nTotal number of distinct parallelograms is {len(sorted_solutions)}.")

solve()
<<<9>>>