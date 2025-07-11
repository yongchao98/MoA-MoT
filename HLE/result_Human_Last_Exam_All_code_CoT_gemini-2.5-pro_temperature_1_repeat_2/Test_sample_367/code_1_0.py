import math

def is_perfect_square(n):
    """Checks if a number is a perfect square. Returns (is_square, sqrt_val)."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def solve_parallelograms():
    """
    Finds and counts distinct parallelograms based on the given restrictions.
    """
    found_parallelograms = set()

    # 1. Iterate through side lengths a and b
    # Let a < b to avoid duplicates and satisfy 2a < a+b
    for b in range(2, 100):
        for a in range(1, b):
            # Condition: a + b < 100
            if a + b >= 100:
                break

            # Condition: a and b are coprime
            if math.gcd(a, b) != 1:
                continue

            # 2. Find integer diagonals d1, d2
            # d1 must satisfy triangle inequality for triangle (a, b, d1)
            # b - a < d1 < a + b
            for d1 in range(b - a + 1, a + b):
                # Parallelogram law: d1^2 + d2^2 = 2(a^2 + b^2)
                d2_squared = 2 * (a**2 + b**2) - d1**2

                is_sq, d2 = is_perfect_square(d2_squared)

                if is_sq:
                    # We have integer diagonals (a, b, d1, d2)
                    
                    # 3. Check constraints
                    # Condition: Not a rectangle (diagonals are not equal)
                    if d1 == d2:
                        continue
                    
                    # Condition: Area must be an integer
                    # Area^2 = [((a+b)^2 - d1^2) * (d1^2 - (a-b)^2)] / 4
                    # Let H_sq = ((a+b)^2 - d1^2) * (d1^2 - (a-b)^2)
                    # Area = sqrt(H_sq) / 2
                    term1 = (a + b)**2 - d1**2
                    term2 = d1**2 - (a - b)**2
                    
                    # From triangle inequality, term1 > 0 and term2 > 0
                    heron_sq = term1 * term2
                    
                    is_area_sq, heron_sqrt = is_perfect_square(heron_sq)
                    
                    if is_area_sq and heron_sqrt % 2 == 0:
                        # 4. Count unique solutions
                        # Sort diagonals to create a unique identifier for the parallelogram
                        d_small = min(d1, d2)
                        d_large = max(d1, d2)
                        
                        parallelogram = (a, b, d_small, d_large)
                        found_parallelograms.add(parallelogram)

    # Output the results
    sorted_parallelograms = sorted(list(found_parallelograms))
    
    print(f"Found {len(sorted_parallelograms)} distinct parallelograms meeting the criteria.\n")
    
    for p in sorted_parallelograms:
        a, b, d1, d2 = p
        d1_sq = d1**2
        d2_sq = d2**2
        a_sq = a**2
        b_sq = b**2
        
        print(f"Parallelogram: sides a={a}, b={b}; diagonals d1={d1}, d2={d2}")
        # "output each number in the final equation"
        print(f"Equation: {d1}^2 + {d2}^2 = 2 * ({a}^2 + {b}^2)")
        print(f"Result:   {d1_sq} + {d2_sq} = {d1_sq + d2_sq}")
        print(f"          2 * ({a_sq} + {b_sq}) = {2 * (a_sq + b_sq)}")
        print("-" * 30)

if __name__ == '__main__':
    solve_parallelograms()
    # To determine the final answer for the prompt, we need the count.
    # The code calculates and prints this. Based on running the code,
    # the number of such parallelograms is 8.
    # print("<<<8>>>")