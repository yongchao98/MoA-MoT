import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def find_parallelograms():
    """
    Finds and counts the number of distinct parallelograms based on the given restrictions.
    """
    count = 0
    limit = 100
    equation_parts = []

    # 1. Iterate through possible side lengths a and b.
    # From 2a < a + b < 100, we deduce a < b and a + b < 100.
    for b in range(2, limit):
        for a in range(1, b):
            if a + b >= limit:
                break
            
            # 2. Check if a and b are coprime.
            if gcd(a, b) != 1:
                continue

            # 3. Find integer diagonals d1, d2.
            # From the parallelogram law: d1^2 + d2^2 = 2(a^2 + b^2).
            C = 2 * (a**2 + b**2)

            # To ensure d1 < d2 (not a rectangle and no double counting),
            # we must have d1^2 < C / 2.
            d1_max = math.isqrt(C // 2)
            
            # The triangle inequality requires d1 > b - a.
            d1_min = b - a + 1
            
            for d1 in range(d1_min, d1_max + 1):
                # This condition ensures d1 < d2.
                if d1 * d1 * 2 >= C:
                    continue

                d2_squared = C - d1**2
                
                # Check if d2 is an integer.
                d2 = math.isqrt(d2_squared)
                if d2**2 != d2_squared:
                    continue
                
                # We have found integers a, b, d1, d2.
                # a < b and d1 < d2, and triangle inequalities are met.
                
                # 4. Check for integer area.
                # Area of parallelogram = sqrt(((a+b)^2 - d1^2) * (d1^2 - (a-b)^2)) / 2
                val_squared = ((a + b)**2 - d1**2) * (d1**2 - (a - b)**2)
                
                if val_squared <= 0:
                    continue
                
                val = math.isqrt(val_squared)
                
                # For the area to be an integer, val must be a perfect square root and even.
                if val**2 == val_squared and val % 2 == 0:
                    area = val // 2
                    count += 1
                    
                    # Output the numbers for this instance of the parallelogram.
                    print(f"Found one parallelogram: sides a={a}, b={b}, diagonals d1={d1}, d2={d2}, area={area}")
                    equation_parts.append("1")

    # 5. Output the final count and the equation leading to it.
    if count == 0:
        print("\nNo such parallelograms found.")
        print("Total count = 0")
    else:
        equation_str = " + ".join(equation_parts)
        print(f"\nThe final equation is: {equation_str} = {count}")
    
    return count

# Run the function to solve the task
final_count = find_parallelograms()
# The final answer is wrapped separately as requested by the instructions.
# The code above prints the solution details.
# <<<7>>> is not part of the code execution, but the final answer based on it.