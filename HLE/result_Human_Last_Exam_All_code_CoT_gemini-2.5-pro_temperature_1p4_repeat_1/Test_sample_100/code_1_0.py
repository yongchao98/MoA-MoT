import math

def solve_integral():
    """
    This function calculates the value of the given integral based on the
    symbolic derivation.
    The integral simplifies to the expression:
    (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
    """
    pi = math.pi
    
    # The final equation is of the form:
    # (n1/d1)*pi^8 + (n2/d2)*pi^2 - (n3/d3)*pi + n4
    
    n1 = 8
    d1 = 15
    n2 = 1
    d2 = 3
    n3 = 1
    d3 = 2
    n4 = 1

    # Calculate the value of the expression
    value = (n1/d1) * pi**8 + (n2/d2) * pi**2 - (n3/d3) * pi + n4

    # Output the symbolic form of the final equation, showing each number
    print("The final simplified expression for the integral is:")
    print("({}/{}) * pi^8 + ({}/{}) * pi^2 - ({}/{}) * pi + {}".format(n1, d1, n2, d2, n3, d3, n4))
    
    # Output the numerical value
    print("\nThe numerical value of the integral is:")
    print(value)

solve_integral()