import math

def solve_parallelograms():
    """
    Finds and counts the number of distinct parallelograms based on a given set of rules.

    The rules for the parallelogram are:
    1. It is not a rhombus (a != b) and not a rectangle.
    2. Side lengths a and b are coprime integers with 2a < a + b < 100.
       This simplifies to a < b and a + b < 100.
    3. The area is an integer.
    4. Both diagonals (d1, d2) have integer lengths.

    The solution is based on these mathematical properties:
    - The Parallelogram Law: d1^2 + d2^2 = 2 * (a^2 + b^2).
    - Area formula: 16 * Area^2 = (4ab)^2 - (d2^2 - d1^2)^2.
    """
    
    count = 0
    
    # Iterate through all possible side lengths a and b.
    # From a < b and a + b < 100, b can go up to 98 (for a=1).
    for b in range(2, 99):
        for a in range(1, b):
            # Condition: a + b < 100
            if a + b >= 100:
                break
            
            # Condition: a and b are coprime
            if math.gcd(a, b) != 1:
                continue

            # From the Parallelogram Law, calculate the required sum of squares for the diagonals.
            s = 2 * (a**2 + b**2)
            
            # Search for integer diagonals d1, d2.
            # To avoid duplicates and unnecessary checks, we iterate d1 up to sqrt(s/2)
            # which ensures d1 <= d2.
            limit_d1 = int(math.sqrt(s / 2))
            for d1 in range(1, limit_d1 + 1):
                d2_squared = s - d1**2
                
                # Check if d2 is an integer
                d2 = int(math.sqrt(d2_squared))
                if d2**2 != d2_squared:
                    continue

                # Now we have a candidate set (a, b, d1, d2).
                
                # Condition: Not a rectangle (d1 != d2).
                # The loop condition d1 <= limit_d1 < d2 ensures d1 != d2.
                if d1 == d2:
                    continue
                
                # Condition: Triangle inequalities must hold for a valid, non-degenerate parallelogram.
                # The diagonals must be greater than the difference of sides and less than the sum.
                if not (b - a < d1 < a + b and b - a < d2 < a + b):
                    continue

                # Condition: The area must be an integer.
                # This is checked by verifying that (4ab)^2 - (d2^2 - d1^2)^2 is a perfect square K^2,
                # and K must be divisible by 4 (since K = 4 * Area).
                
                d_diff_sq = d2_squared - d1**2
                area_term_sq = (4 * a * b)**2 - d_diff_sq**2
                
                if area_term_sq < 0:
                    continue
                
                area_term_root = int(math.sqrt(area_term_sq))
                
                if area_term_root**2 == area_term_sq and area_term_root % 4 == 0:
                    # All conditions are satisfied. We have found a distinct parallelogram.
                    count += 1
                    area = area_term_root // 4
                    
                    # Output the details of the found parallelogram and the verification equation.
                    print(f"Found parallelogram: a={a}, b={b}, d1={d1}, d2={d2}, Area={area}")
                    print(f"  Equation check: {d1}^2 + {d2}^2 = {d1**2 + d2**2}")
                    print(f"  Equation check: 2*(a^2 + b^2) = 2*({a**2} + {b**2}) = {s}")


    print(f"\nTotal number of distinct parallelograms is: {count}")

solve_parallelograms()
<<<9>>>