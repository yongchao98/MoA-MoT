import math

def check_if_composite(n):
    """
    This function implements a semi-decision procedure for the set of composite numbers.
    It searches for factors x and y to satisfy the existential formula.
    If the number n is composite, it finds the factors and halts, returning True.
    If n is not composite (i.e., prime or <= 1), the search will complete and it will return False.
    This corresponds to the set being recursive (decidable), and thus also RE.
    
    The defining formula is:
    exists x,y,w_x,w_y ( P(x) & P(y) & P(w_x) & P(w_y) & x=w_x+2 & y=w_y+2 & v=x*y )
    which simplifies to: exists x, y in N, (x>=2 & y>=2 & v=x*y)
    
    Here, v is the input number 'n'.
    """
    if not isinstance(n, int) or n <= 1:
        return False

    # The search for x and y corresponds to resolving the existential quantifiers.
    # We can bound the search for x up to sqrt(n).
    for x in range(2, int(math.sqrt(n)) + 1):
        if n % x == 0:
            y = n // x
            # We found x and y such that n = x * y, and x,y >= 2
            # Let's print the numbers in the final satisfied equation.
            print(f"Checking n = {n}: Found solution to n = x * y. Equation: {n} = {x} * {y}")
            return True
            
    print(f"Checking n = {n}: No solution found in the search space.")
    return False

# Demonstrate for a few numbers.
numbers_to_test = [4, 6, 7, 9, 13, 100]
for num in numbers_to_test:
    is_composite = check_if_composite(num)
    print(f"Is {num} composite? {is_composite}\n")
