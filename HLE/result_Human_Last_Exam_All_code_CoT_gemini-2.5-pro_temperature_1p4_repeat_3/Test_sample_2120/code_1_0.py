import sys

def solve():
    """
    This function calculates the average value of the complex coordinates z
    where the inverse of the gradient of S(z) approaches zero.
    """
    # As outlined in the plan, the problem reduces to finding the average of the roots
    # of the polynomial in the denominator of the B(z) equation.
    # The polynomial is: 4*z^4 - z^3 + z^2 + 1 = 0.
    
    # Coefficients of the polynomial a*z^4 + b*z^3 + c*z^2 + d*z + e = 0
    a = 4  # coefficient of z^4
    b = -1 # coefficient of z^3
    # c = 1, d = 0, e = 1 are not needed for the sum of roots.
    
    # The degree of the polynomial
    degree = 4
    
    # According to Vieta's formulas, the sum of the roots of a polynomial is -b/a.
    sum_of_roots = -b / a
    
    # The average value of the roots is their sum divided by the number of roots (the degree).
    average_of_roots = sum_of_roots / degree
    
    print("The critical coordinates 'z' are the roots of the polynomial:")
    print(f"{a}*z^4 + ({b})*z^3 + 1*z^2 + 1 = 0")
    print("")
    print("To find the average of these roots, we use Vieta's formulas.")
    print(f"The degree of the polynomial is {degree}.")
    print(f"The coefficient of z^4 (a_4) is {a}.")
    print(f"The coefficient of z^3 (a_3) is {b}.")
    print("")
    
    # Outputting each number in the final equation as requested
    print(f"The sum of the roots is calculated as -a_3 / a_4:")
    print(f"Sum = -({b}) / {a} = {sum_of_roots}")
    
    print(f"The average of the roots is the sum divided by the degree:")
    print(f"Average = {sum_of_roots} / {degree} = {average_of_roots}")

solve()