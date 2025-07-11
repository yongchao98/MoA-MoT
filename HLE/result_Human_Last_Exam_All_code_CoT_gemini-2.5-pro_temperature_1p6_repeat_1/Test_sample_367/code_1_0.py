import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve_parallelograms():
    """
    Finds the number of distinct parallelograms based on the given restrictions.
    """
    solutions = set()
    
    # 1. Generate unique Primitive Pythagorean Triples (PPTs)
    # The hypotenuse z must divide a*b. Since a+b < 100, the max a*b is < 49*50=2450.
    # z = m^2 + n^2, so m is bounded by sqrt(2450) ~= 49.5.
    unique_ppts = set()
    for m in range(2, 50):
        for n in range(1, m):
            if (m - n) % 2 == 1 and math.gcd(m, n) == 1:
                x = m*m - n*n
                y = 2*m*n
                z = m*m + n*n
                # Store a canonical form (min, max, hyp) to keep the set unique
                unique_ppts.add((min(x, y), max(x, y), z))

    # We need to test both orientations for the angle, i.e., cos(theta) = x/z and cos(theta) = y/z
    all_oriented_ppts = []
    for x_s, y_s, z in sorted(list(unique_ppts)):
        all_oriented_ppts.append((x_s, y_s, z))
        all_oriented_ppts.append((y_s, x_s, z))

    # 2. Iterate through PPTs and find qualifying (a, b) pairs
    for x_prim, y_prim, z_prim in all_oriented_ppts:
        
        # ab = M * z_prim. The max possible value for ab is just under 2500.
        max_ab = 49 * 50
        for M in range(1, (max_ab // z_prim) + 1):
            target_product = M * z_prim

            # 3. Find factor pairs (a,b) of the target product
            for a in range(1, int(math.sqrt(target_product)) + 1):
                if target_product % a == 0:
                    b = target_product // a

                    # 4. Apply all problem restrictions on sides a and b
                    if a >= b:  # Condition 1: a != b (enforced by a < b)
                        continue
                    if a + b >= 100:  # Condition 2: a+b < 100
                        continue
                    if math.gcd(a, b) != 1:  # Condition 2: a and b are coprime
                        continue

                    # 5. Check if both diagonals are integers
                    # d^2 = a^2+b^2 +/- 2*a*b*cos(theta) = a^2+b^2 +/- 2*(M*z_prim)*(x_prim/z_prim)
                    term = 2 * M * x_prim
                    A_val = a*a + b*b
                    
                    d1_sq = A_val - term
                    d2_sq = A_val + term
                    
                    # Also, the parallelogram must not be a rectangle (diagonals are equal).
                    # This happens if term=0, which requires x_prim=0, not possible for a PPT.
                    if is_perfect_square(d1_sq) and is_perfect_square(d2_sq):
                        d1 = int(math.sqrt(d1_sq))
                        d2 = int(math.sqrt(d2_sq))
                        
                        # Store a canonical representation of the parallelogram's geometry
                        solution_tuple = (a, b, min(d1, d2), max(d1, d2))
                        solutions.add(solution_tuple)

    print(len(solutions))

solve_parallelograms()