import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = math.isqrt(n)
    return x * x == n

def solve():
    """
    Finds the number of distinct parallelograms satisfying the given conditions.
    """
    found_parallelograms = set()
    
    # Loop through possible side lengths b and a
    # Conditions: a < b and a + b < 100
    for b in range(2, 99):
        # a < b and a < 100 - b
        for a in range(1, min(b, 100 - b)):
            
            # Condition 2: a and b are coprime
            if math.gcd(a, b) != 1:
                continue

            # Condition 1: Not a rectangle
            # If a^2 + b^2 is a perfect square, it's a rectangle.
            if is_perfect_square(a * a + b * b):
                continue
            
            # Pre-calculate terms for efficiency
            ab_sq = (a * b) ** 2
            sum_of_sq_sides = 2 * (a * a + b * b)

            # Iterate through possible lengths for the first diagonal, d1
            # Triangle inequality for triangle (a, b, d1): d1 < a+b and d1 > b-a
            for d1 in range(b - a + 1, a + b):
                
                # From d1^2 + d2^2 = 2*(a^2 + b^2), calculate d2^2
                d2_sq = sum_of_sq_sides - d1 * d1
                
                # d2 must be an integer, so d2_sq must be a perfect square
                if d2_sq <= 0 or not is_perfect_square(d2_sq):
                    continue
                
                d2 = math.isqrt(d2_sq)

                # To count unique parallelograms, we enforce an order d1 < d2
                # d1=d2 implies a rectangle, which is already excluded.
                if d1 >= d2:
                    continue

                # Condition 3: Area must be an integer
                # This is equivalent to (a*b)^2 - q^2 being a perfect square,
                # where q = (d2^2 - d1^2) / 4.
                
                q_numerator = d2 * d2 - d1 * d1
                
                # d1 and d2 must have the same parity for q to be an integer.
                # The sum of squares identity guarantees this.
                # So d2^2 - d1^2 is always divisible by 4.
                if q_numerator % 4 != 0:
                    continue # Should not happen, but a good safeguard.
                
                q = q_numerator // 4
                
                area_sq = ab_sq - q * q
                
                if area_sq > 0 and is_perfect_square(area_sq):
                    # All conditions are met. We have found a valid parallelogram.
                    # The tuple (a, b, d1, d2) uniquely defines it.
                    parallelogram = (a, b, d1, d2)
                    found_parallelograms.add(parallelogram)

    print(len(found_parallelograms))

solve()