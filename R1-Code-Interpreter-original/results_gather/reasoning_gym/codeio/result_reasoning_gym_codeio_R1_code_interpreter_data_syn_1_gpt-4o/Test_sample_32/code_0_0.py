# import necessary packages
from sympy import symbols, Eq, limit, solve

# main function
def main_solution(delta: float):
    # Define the symbols
    x, k = symbols('x k')
    
    # Define the function f(x)
    def f(a):
        if a <= 5:
            return k * x + 1
        elif a > 5:
            return 3 * x - 5
    
    # Calculate the left and right hand limits at x = 5
    LHL = limit(f(5 - delta), x, 5, '-')
    RHL = limit(f(5 + delta), x, 5, '+')
    
    # Create the equation for continuity
    eq1 = Eq(LHL, RHL)
    
    # Solve for k
    solution = solve(eq1, k)
    
    # Return the value of k as a float
    return float(solution[0])

# Test the function with a small delta
delta = 0.01
k_value = main_solution(delta)
print(k_value)