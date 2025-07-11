import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve():
    """
    Finds the number of distinct parallelograms satisfying the given restrictions.
    """
    found_parallelograms = set()

    # 1. Iterate through all possible side lengths a and b.
    # The condition 2a < a + b < 100 implies a < b and a + b < 100.
    # The largest possible 'a' is 48 (e.g., a=48, b=51 -> a+b=99).
    for a in range(1, 50):
        for b in range(a + 1, 100 - a):
            
            # 2. Check for coprime sides.
            if gcd(a, b) != 1:
                continue

            # This is the required sum of squares of the diagonals.
            diag_sq_sum = 2 * (a**2 + b**2)

            # 3. Iterate through possible integer diagonals d1.
            # d1 must satisfy the triangle inequality: b-a < d1 < a+b.
            for d1 in range(b - a + 1, a + b):
                
                # 4a. Check if it's a rectangle (d1^2 = a^2 + b^2).
                if d1**2 == a**2 + b**2:
                    continue

                # 4b. Check if the second diagonal d2 is an integer.
                d2_sq = diag_sq_sum - d1**2
                if d2_sq <= 0 or not is_perfect_square(d2_sq):
                    continue
                
                d2 = int(math.sqrt(d2_sq))

                # d2 must also satisfy the triangle inequality.
                if not (b - a < d2 < a + b):
                    continue

                # 4c. Check if the area is an integer.
                # This is equivalent to checking if ((a+b)^2 - d1^2) * (d1^2 - (b-a)^2) is a perfect square.
                area_term_1 = (a + b)**2 - d1**2
                area_term_2 = d1**2 - (b - a)**2
                
                if is_perfect_square(area_term_1 * area_term_2):
                    # A valid parallelogram is found.
                    # Store it in a canonical form to count unique ones.
                    canonical_form = (a, b, min(d1, d2), max(d1, d2))
                    found_parallelograms.add(canonical_form)

    # 5. Print the results.
    if not found_parallelograms:
        print("No parallelograms were found that satisfy all the given conditions.")
    else:
        print(f"Found {len(found_parallelograms)} distinct parallelograms satisfying the conditions.")
        print("For each parallelogram (a, b) with diagonals (d1, d2), the parallelogram law holds:")
        
        sorted_grams = sorted(list(found_parallelograms))
        for a_sol, b_sol, d1_sol, d2_sol in sorted_grams:
            lhs = d1_sol**2 + d2_sol**2
            rhs = 2 * (a_sol**2 + b_sol**2)
            print(f"Sides: ({a_sol}, {b_sol}), Diagonals: ({d1_sol}, {d2_sol})")
            print(f"Equation: {d1_sol}**2 + {d2_sol}**2 = {lhs} = 2*({a_sol}**2 + {b_sol}**2)")
            print("-" * 20)
            
    print(f"\nTotal number of distinct parallelograms: {len(found_parallelograms)}")

if __name__ == '__main__':
    solve()