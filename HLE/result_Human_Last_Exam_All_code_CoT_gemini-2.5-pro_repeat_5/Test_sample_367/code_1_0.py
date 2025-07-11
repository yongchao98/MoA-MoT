import math

def is_perfect_square(n):
    """Checks if a number is a perfect square and returns the root if it is."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def solve():
    """
    Finds the number of distinct parallelograms satisfying the given restrictions.
    """
    count = 0
    found_parallelograms = []

    # Iterate through possible side lengths a and b
    # Assume a < b to count each pair of side lengths only once.
    for b in range(2, 100):
        for a in range(1, b):
            # Condition 2: 2a < a + b < 100  (a < b implies 2a < a+b)
            if a + b >= 100:
                break
            
            # Condition 2: a and b are coprime
            if gcd(a, b) != 1:
                continue
            
            # From the parallelogram law: d1^2 + d2^2 = 2(a^2 + b^2)
            C = 2 * (a*a + b*b)
            
            # Find integer solutions (d1, d2) for d1^2 + d2^2 = C
            # Loop d1 up to sqrt(C/2) to ensure d1 <= d2
            limit = int(math.sqrt(C / 2))
            for d1 in range(1, limit + 1):
                d2_sq = C - d1*d1
                is_sq, d2 = is_perfect_square(d2_sq)
                
                if is_sq:
                    # Found a potential pair of diagonals (d1, d2)
                    
                    # Condition 1: Not a rectangle (d1 != d2)
                    if d1 == d2:
                        continue
                    
                    # Condition for non-degenerate triangles (a,b,d1) and (a,b,d2)
                    # This means a+b > d2 > d1 > b-a
                    if not (a + b > d2 and d1 > b - a):
                        continue

                    # Condition 3: Area is an integer.
                    # This holds if the perimeters a+b+d1 and a+b+d2 are even.
                    # This implies d1 and d2 must have the same parity as a+b.
                    ab_sum_is_even = ((a + b) % 2 == 0)
                    d1_is_even = (d1 % 2 == 0)
                    d2_is_even = (d2 % 2 == 0)

                    # Case 1: a and b have different parity, so a+b is odd.
                    # Diagonals d1 and d2 must both be odd.
                    if not ab_sum_is_even and (d1_is_even or d2_is_even):
                        continue
                    
                    # Case 2: a and b are both odd (coprime means not both even), so a+b is even.
                    # Diagonals d1 and d2 must both be even.
                    if ab_sum_is_even and (not d1_is_even or not d2_is_even):
                        continue
                    
                    # All conditions are met for the parallelogram (a, b, d1, d2).
                    count += 1
                    # This problem has an equation like structure where each found set is a solution
                    print(f"Solution {count}: a={a}, b={b}, diagonals=({d1}, {d2})")
    
    print(f"\nThe total number of distinct parallelograms is: {count}")

solve()