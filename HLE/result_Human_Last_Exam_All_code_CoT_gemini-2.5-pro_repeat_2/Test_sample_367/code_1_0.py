import math

def solve():
    """
    Calculates the number of distinct parallelograms satisfying the given criteria.
    """
    
    def is_perfect_square(n):
        """Checks if a number is a perfect square. Returns (bool, sqrt_val)."""
        if n < 0:
            return False, -1
        if n == 0:
            return True, 0
        x = int(math.sqrt(n))
        return x * x == n, x

    # Use a set to store unique solutions to avoid double counting.
    # A distinct parallelogram is defined by its sides (a,b) and one diagonal d1.
    solutions = set()

    # 1. Iterate through possible side lengths a and b.
    # The condition a+b < 100 and a < b implies 2*a + 1 < 100, so a <= 49.
    for a in range(1, 50):
        # The upper bound for b is 99-a from a+b < 100.
        for b in range(a + 1, 100 - a):
            
            # Condition: a and b must be coprime.
            if math.gcd(a, b) != 1:
                continue

            # 2. Find possible diagonals using the parallelogram law: d1^2 + d2^2 = 2(a^2 + b^2).
            # To avoid duplicates and rectangles, we enforce d1 < d2.
            # This implies d1^2 < a^2 + b^2.
            # The triangle inequality for sides (a, b, d1) requires d1 > b - a.
            
            N = 2 * (a*a + b*b)
            limit_d1_sq = a*a + b*b
            
            # We search for d1 in the valid range.
            d1_start = b - a + 1
            d1_end = int(math.sqrt(limit_d1_sq))

            for d1 in range(d1_start, d1_end + 1):
                # Ensure d1^2 < a^2+b^2 to exclude rectangles.
                if d1 * d1 >= limit_d1_sq:
                    continue

                d2_sq = N - d1*d1
                is_sq, d2 = is_perfect_square(d2_sq)

                if is_sq:
                    # We have a candidate parallelogram with sides (a,b) and diagonals (d1,d2).
                    # 3. Check for the integer area condition.
                    # Area = (1/2) * sqrt(d1^2*d2^2 - (a^2-b^2)^2).
                    
                    term1 = d1*d1 * d2_sq
                    term2 = (a*a - b*b)**2
                    
                    if term1 <= term2:
                        continue
                        
                    area_sq_num = term1 - term2
                    is_sq_area, area_num = is_perfect_square(area_sq_num)
                    
                    if is_sq_area:
                        # For Area to be an integer, area_num must be even.
                        if area_num % 2 == 0:
                            # All conditions are met. This is a valid parallelogram.
                            # We add (a,b,d1) to our set of solutions.
                            solution_tuple = (a, b, d1, d2, area_num // 2)
                            solutions.add(solution_tuple)

    # Output the details of each found parallelogram
    # This fulfills the "output each number in the final equation" instruction.
    sorted_solutions = sorted(list(solutions))
    for sol in sorted_solutions:
        a, b, d1, d2, area = sol
        print(f"Found parallelogram: a={a}, b={b}, diagonals=({d1},{d2}), area={area}")

    # The final answer is the total number of distinct parallelograms found.
    print(f"\nTotal number of distinct parallelograms:")
    print(len(solutions))
    
    return len(solutions)

# Run the solver and store the final count.
final_count = solve()
# The final answer is presented in the required format.
# print(f"<<<{final_count}>>>")