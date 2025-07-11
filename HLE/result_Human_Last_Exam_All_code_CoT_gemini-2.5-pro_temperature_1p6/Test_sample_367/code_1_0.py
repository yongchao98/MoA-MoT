import math

def main():
    """
    This function calculates the number of distinct parallelograms satisfying the given criteria.
    A distinct parallelogram is defined by a unique set of side lengths {a, b} and diagonal lengths {d1, d2}.
    """
    count = 0
    # 1. Iterate through possible side lengths a and b.
    # To satisfy 2a < a+b, we must have a < b. We loop b from 2 up to 98.
    for b in range(2, 100):
        # a is the shorter side.
        for a in range(1, b):
            # To satisfy a+b < 100. If violated, inner loop can stop.
            if a + b >= 100:
                break
            
            # Sides must be coprime.
            if math.gcd(a, b) != 1:
                continue

            # Parallelogram is not a rhombus (a != b), guaranteed by loop ranges.

            # 2. Use the parallelogram law: d1^2 + d2^2 = 2*(a^2 + b^2).
            c_val = 2 * (a**2 + b**2)

            # 3. Iterate through possible diagonal lengths d1.
            # Enforce d1 < d2 (for non-rectangles and unique counting), which implies d1^2 < c_val / 2.
            # 4. Apply triangle inequality: d1 > b - a.
            d1_lower_bound = b - a + 1
            d1_upper_bound = int(math.sqrt(c_val / 2.0))

            for d1 in range(d1_lower_bound, d1_upper_bound + 1):
                d2_squared = c_val - d1**2
                
                # Check if d2 is an integer.
                d2 = int(math.sqrt(d2_squared))
                if d2**2 == d2_squared:
                    # Now we have integer diagonals d1, d2.
                    # Loop structure ensures d1 < d2, so it's not a rectangle.
                    if d1 == d2:
                        continue
                        
                    # Apply second part of triangle inequality: d2 < a + b.
                    if d2 < a + b:
                        # 5. Area is guaranteed to be an integer.
                        # All conditions met. We found a valid parallelogram.
                        count += 1
                        
    print(count)

if __name__ == '__main__':
    main()