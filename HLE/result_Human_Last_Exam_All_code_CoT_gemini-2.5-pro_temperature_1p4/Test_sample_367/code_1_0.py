import math

def is_perfect_square(n):
    """
    Checks if a non-negative integer n is a perfect square.
    Returns a tuple (is_square, sqrt_of_n).
    """
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def find_parallelograms():
    """
    Finds and counts distinct parallelograms based on the given restrictions.
    A distinct parallelogram is uniquely identified by its side lengths {a, b}
    and diagonal lengths {d1, d2}.
    """
    solutions = []
    
    # Condition: 2a < a + b < 100, which means a < b and a + b < 100.
    # The largest possible 'a' is 48 (with b=51, a+b=99 fails gcd, try b=50 fails gcd, b=49 fails a<b).
    # a=48, b=49 -> a+b=97. So a goes up to 48. range(1,49).
    for a in range(1, 50):
        # b must be > a and a + b < 100
        for b in range(a + 1, 100 - a):
            
            # Condition: a and b must be coprime.
            if math.gcd(a, b) != 1:
                continue

            # Based on a geometric construction, we search for an integer x.
            # (x, y, b) must form a Pythagorean triple.
            # We iterate through possible integer values for x, where 1 <= x < b.
            # x=0 would imply a rectangle, which is disallowed.
            # x=b would imply a degenerate parallelogram with zero area.
            for x in range(1, b):
                # Calculate y^2 from b^2 = x^2 + y^2
                y_squared = b*b - x*x
                is_sq_y, y = is_perfect_square(y_squared)

                if is_sq_y:
                    # Found a valid Pythagorean triple (x, y, b).
                    # Now, check if the diagonals are integers.
                    # d1^2 = (a-x)^2 + y^2 = a^2 - 2ax + b^2
                    # d2^2 = (a+x)^2 + y^2 = a^2 + 2ax + b^2
                    
                    d1_squared = a*a + b*b - 2*a*x
                    is_sq_d1, d1 = is_perfect_square(d1_squared)
                    
                    if is_sq_d1:
                        d2_squared = a*a + b*b + 2*a*x
                        is_sq_d2, d2 = is_perfect_square(d2_squared)
                        
                        if is_sq_d2:
                            # All conditions are met. A valid parallelogram is found.
                            # It's not a rhombus (a < b) and not a rectangle (x > 0 implies d1 != d2).
                            # We store the solution, ensuring d1 < d2 for uniqueness.
                            if d1 < d2:
                                solutions.append((a, b, d1, d2))
                            else:
                                solutions.append((a, b, d2, d1))
    
    # Use a set to count only the unique parallelograms (defined by their dimensions).
    unique_solutions = sorted(list(set(solutions)))
    
    # Print the results
    count = len(unique_solutions)
    print(f"Found {count} distinct parallelograms that satisfy all the conditions.")
    
    if count > 0:
        print("\nHere are the details for each parallelogram:")
        for a, b, d1, d2 in unique_solutions:
            print(f"\n- Parallelogram with sides a={a}, b={b} and diagonals d1={d1}, d2={d2}")
            # Output the parallelogram law equation with the numbers plugged in
            print(f"  The final equation is: 2 * ({a}**2 + {b}**2) = {d1}**2 + {d2}**2")
            print(f"  Calculation: 2 * ({a*a} + {b*b}) = {d1*d1} + {d2*d2}")
            print(f"  Result: {2 * (a*a + b*b)} = {d1*d1 + d2*d2}")
            
    # Return the final count as requested.
    return count

if __name__ == '__main__':
    final_answer = find_parallelograms()
    print(f"\n<<<Number of distinct parallelograms: {final_answer}>>>")