def check_composite_definition(n):
    """
    Checks if a number n is composite based on its Diophantine representation.
    A number n is composite if there exist natural numbers x and y such that n = (x+2)*(y+2).
    This polynomial equation Q(n, x, y) = n - (x+2)*(y+2) = 0 has a solution in N
    if and only if n is composite. (This covers all composites >= 4).
    We will search for non-negative integers x and y that satisfy the equation.
    
    The corresponding L-formula is: exists x, y (P(x) AND P(y) AND n - (x+2)*(y+2) = 0)
    assuming P defines non-negative integers.
    """
    if not isinstance(n, int) or n < 4:
        return False, None, None

    # We need to find if there exists a solution (x, y) in natural numbers.
    # We can bound the search. Since x, y >= 0, we have x+2 >= 2 and y+2 >= 2.
    # So, n = (x+2)*(y+2) >= 4.
    # We can assume x+2 <= sqrt(n), so x <= sqrt(n) - 2.
    limit = int(n**0.5)
    for x_plus_2 in range(2, limit + 1):
        if n % x_plus_2 == 0:
            y_plus_2 = n // x_plus_2
            x = x_plus_2 - 2
            y = y_plus_2 - 2
            # We found a solution in non-negative integers
            return True, x, y
            
    return False, None, None

# Test cases
test_numbers = [4, 5, 6, 9, 13, 25, 33]

print("This program demonstrates that the set of composite numbers, a recursively enumerable set,")
print("is existentially definable. A number n is composite iff the equation")
print("n = (x+2)*(y+2) has a solution for some natural numbers x and y.")
print("-" * 20)

for n in test_numbers:
    is_comp, x_sol, y_sol = check_composite_definition(n)
    if is_comp:
        print(f"For n = {n}: It IS composite.")
        print(f"The equation {n} = (x+2)*(y+2) has a solution (x, y) = ({x_sol}, {y_sol}).")
        # Final equation output for the solution
        print(f"Final equation: {n} - ({x_sol} + 2) * ({y_sol} + 2) = {n - (x_sol + 2) * (y_sol + 2)}")
    else:
        print(f"For n = {n}: It is NOT composite.")
        print(f"The equation {n} = (x+2)*(y+2) has no solution for natural numbers x, y.")
    print("-" * 20)