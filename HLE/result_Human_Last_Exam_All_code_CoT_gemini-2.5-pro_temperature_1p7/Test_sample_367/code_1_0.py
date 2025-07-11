import math

def solve():
    """
    This function calculates the number of distinct parallelograms satisfying the given criteria.
    1.  The parallelogram is neither a rhombus (a!=b) nor a rectangle (d1!=d2).
    2.  Side lengths a and b are coprime integers with 2a < a + b < 100.
    3.  The area of the parallelogram is an integer.
    4.  The lengths of both diagonals (d1, d2) are integers.
    """
    
    def is_positive_perfect_square(n):
        if n <= 0:
            return False
        x = int(math.sqrt(n))
        return x * x == n

    count = 0
    first_solution_details = None

    # Loop for side length a
    for a in range(1, 50):
        # Loop for side length b
        # a < b is derived from 2a < a + b
        for b in range(a + 1, 100):
            # Constraint: a + b < 100
            if a + b >= 100:
                break

            # Constraint: a and b are coprime
            if math.gcd(a, b) != 1:
                continue
            
            # From the Parallelogram Law: d1^2 + d2^2 = 2(a^2 + b^2)
            C = 2 * (a**2 + b**2)
            
            # Iterate d1 up to sqrt(C/2). This ensures d1 <= d2.
            d1_limit = int(math.sqrt(C / 2.0))
            for d1 in range(1, d1_limit + 1):
                d2_squared = C - d1**2
                
                if is_positive_perfect_square(d2_squared):
                    d2 = int(math.sqrt(d2_squared))

                    # Constraint: Not a rectangle (d1 != d2).
                    if d1 == d2:
                        continue

                    # A non-degenerate parallelogram requires d1 > b-a
                    if not (d1 > b - a):
                        continue
                    
                    # From Area^2 = (a*b)^2 - K^2, where K = (d2^2 - d1^2)/4
                    # For Area to be an integer, (a*b)^2 - K^2 must be a perfect square.
                    K_numerator = d2**2 - d1**2
                    
                    # d1 and d2 have the same parity, so d1^2-d2^2 is a multiple of 4.
                    K = K_numerator // 4
                    
                    area_squared = (a * b)**2 - K**2
                    
                    if is_positive_perfect_square(area_squared):
                        count += 1
                        if first_solution_details is None:
                             area = int(math.sqrt(area_squared))
                             first_solution_details = (a, b, d1, d2, K, area)
                             
    print(f"Total number of distinct parallelograms is {count}.")
    
    if first_solution_details:
        a, b, d1, d2, K, area = first_solution_details
        print("\n--- Example of a valid parallelogram ---")
        print(f"Sides: a = {a}, b = {b}")
        print(f"Diagonals: d1 = {d1}, d2 = {d2}")
        print("\nVerification of conditions:")
        print("1. Diagonal length equation:")
        print(f"   {d1}^2 + {d2}^2 = {d1**2} + {d2**2} = {d1**2 + d2**2}")
        print(f"   2*({a}^2 + {b}^2) = 2*({a**2} + {b**2}) = {2 * (a**2 + b**2)}")
        
        print("\n2. Integer area equation:")
        print(f"   K = ({d2}^2 - {d1}^2) / 4 = ({d2**2} - {d1**2}) / 4 = {K}")
        print(f"   Area^2 = ({a}*{b})^2 - K^2 = {a*b}^2 - {K}^2 = {(a*b)**2 - K**2}")
        print(f"   Area = sqrt({(a*b)**2 - K**2}) = {area}")
        print("   The area is an integer, so the condition is satisfied.")

solve()
<<<11>>>