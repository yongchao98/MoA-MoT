import math

def get_function_equation(x):
    """
    Returns the string equation for the part of the function f(x)
    that is relevant for a given x >= 1.
    """
    if x < 1:
        # Get the function value at x=1
        k = 0 # x=1 is 2^0
        n = 0
        slope = -2
        y_intercept = 2**(2*n + 2)
        val = slope * 1 + y_intercept + 1
        return f"f(x) = {val} for x < 1"

    # Find which interval [2^k, 2^(k+1)] x is in
    k = math.floor(math.log2(x))
    
    # Add 1 to avoid fixed points in R
    shift = 1 
    
    if k % 2 == 0:
        # k = 2n, Interval is [2^(2n), 2^(2n+1)]
        n = k // 2
        slope = -2
        # f(x) = -2*x + C + 1
        # f(2^(2n+1)) = 0 + 1 => -2*2^(2n+1) + C + 1 = 1 => C = 2*2^(2n+1) = 2^(2n+2)
        y_intercept_term = f"2**(2*{n}+2)"
        return f"f(x) = {slope}*x + {y_intercept_term} + {shift}   (on interval [2**{k}, 2**{k+1}])"
    else:
        # k = 2n+1, Interval is [2^(2n+1), 2^(2n+2)]
        n = (k - 1) // 2
        slope = 4
        # f(x) = 4*x + C + 1
        # f(2^(2n+1)) = 0 + 1 => 4*2^(2n+1) + C + 1 = 1 => C = -4*2^(2n+1) = -2^(2n+3)
        y_intercept_term = f"-2**(2*{n}+3)"
        return f"f(x) = {slope}*x + {y_intercept_term} + {shift}   (on interval [2**{k}, 2**{k+1}])"

def print_function_definition():
    """
    Prints the piecewise definition of the function f(x).
    """
    print("The smallest possible nonzero number of fixed points is 1.")
    print("A function f(x) whose extension has exactly one fixed point in the remainder can be constructed piecewise.")
    print("For x < 1, f(x) is constant, f(x) = f(1).")
    print("For x >= 1, the function is defined on intervals [2**k, 2**(k+1)] for k=0,1,2,...")
    print("\nIf k is even (k=2n):")
    print("On [2**(2n), 2**(2n+1)], the equation is of the form f(x) = a*x + b.")
    a = -2
    b_part1 = 2
    b_part2 = 2
    b_shift = 1
    print(f"The number for 'a' is {a}.")
    print(f"The number for the first part of 'b' is {b_part1}.")
    print(f"The number for the exponent of the first part of 'b' is {b_part2}.")
    print(f"The number for the constant shift is {b_shift}.")
    print("The equation is: f(x) = -2*x + 2**(2*n+2) + 1")

    print("\nIf k is odd (k=2n+1):")
    print("On [2**(2n+1), 2**(2n+2)], the equation is of the form f(x) = a*x + b.")
    a = 4
    b_part1 = -2
    b_part2 = 3
    b_shift = 1
    print(f"The number for 'a' is {a}.")
    print(f"The number for the first part of 'b' is {b_part1}.")
    print(f"The number for the exponent of the first part of 'b' is {b_part2}.")
    print(f"The number for the constant shift is {b_shift}.")
    print("The equation is: f(x) = 4*x - 2**(2*n+3) + 1")

if __name__ == '__main__':
    print_function_definition()
