import math
from collections import defaultdict

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def main():
    """
    Finds the number of distinct parallelograms based on the given restrictions.
    """
    # Step 1: Find all numbers that are sums of two squares in multiple ways.
    # a+b < 100 => max b is 98 (for a=1), max a is 49 (for b=50).
    # Max N = 98^2 + 1^2 = 9605. So we need to check x, y up to sqrt(9605) approx 98.
    # Let's use a slightly larger bound for safety.
    limit = 100
    sums_of_squares = defaultdict(list)
    for x in range(1, limit):
        for y in range(x + 1, limit):
            n = x**2 + y**2
            sums_of_squares[n].append((x, y))

    # Filter for numbers with more than one representation
    multiple_reps = {n: pairs for n, pairs in sums_of_squares.items() if len(pairs) > 1}

    # Step 2: Iterate through the candidates and check all conditions.
    count = 0
    found_parallelograms = []
    
    for n, pairs in multiple_reps.items():
        # Iterate through all combinations of two distinct pairs for a given sum N
        for i in range(len(pairs)):
            for j in range(len(pairs)):
                if i == j:
                    continue

                # Assign one pair as sides (a,b) and the other as (p,q)
                a, b = pairs[i]
                p, q = pairs[j]
                
                # Condition 2: Check side restrictions for (a,b)
                # a < b is guaranteed by loop structure
                # a+b < 100
                # gcd(a,b) == 1
                if (a + b >= 100) or (math.gcd(a, b) != 1):
                    continue

                # Area-related condition from our derivation
                # We need (ab)^2 - (pq)^2 to be a perfect square, which requires ab > pq
                ab = a * b
                pq = p * q
                
                if ab <= pq:
                    continue
                
                area_sq_component = ab**2 - pq**2
                
                # Condition 3 & 4: Area and Diagonals are integers
                if is_perfect_square(area_sq_component):
                    count += 1
                    area = int(math.sqrt(area_sq_component))
                    
                    # Also calculate diagonals for completeness and printing
                    m = pq
                    d1_sq = a**2 + b**2 - 2 * m
                    d2_sq = a**2 + b**2 + 2 * m
                    
                    if is_perfect_square(d1_sq) and is_perfect_square(d2_sq):
                        d1 = int(math.sqrt(d1_sq))
                        d2 = int(math.sqrt(d2_sq))
                        found_parallelograms.append({
                            "sides": (a,b), "others": (p,q), 
                            "area": area, "diagonals": (d1, d2)
                        })

    print(f"Found {count} distinct parallelograms satisfying all conditions.\n")
    if found_parallelograms:
        print("Details of the found parallelograms:")
        for k, p in enumerate(found_parallelograms):
            sides = p["sides"]
            others = p["others"]
            diags = p['diagonals']
            area = p['area']
            
            print(f"  Parallelogram {k+1}:")
            print(f"    - Sides (a,b): {sides}")
            print(f"    - Paired Sum of Squares (p,q): {others}")
            print(f"    - Area: {area}")
            print(f"    - Diagonals: {diags}")

if __name__ == "__main__":
    main()
<<<7>>>