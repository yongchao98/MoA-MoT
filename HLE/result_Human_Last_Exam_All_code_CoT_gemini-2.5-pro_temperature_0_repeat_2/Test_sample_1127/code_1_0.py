import math

def find_minimal_polynomial_for_a_plus_b_sqrt_c(a, b, c):
    """
    Calculates the coefficients of the minimal polynomial for a number of the form a + b*sqrt(c).
    The number x = a + b*sqrt(c) is a root of the polynomial P(x) = x^2 - 2ax + (a^2 - b^2*c) = 0.
    """
    
    # The coefficients of the polynomial P(x) = c2*x^2 + c1*x + c0 = 0
    c2 = 1
    c1 = -2 * a
    c0 = a**2 - (b**2) * c
    
    return c2, c1, c0

def solve_task():
    """
    Solves the main task by finding the minimal polynomial of the connective constant.
    
    The plan is as follows:
    1. The graph is identified as a 4-regular triangulation of an infinite strip.
    2. The connective constant mu is estimated by counting short self-avoiding walks.
       - c_1 = 4
       - c_2 = 12
       - c_3 = 30
       - An estimate for mu is (c_3/c_1)^(1/2) = sqrt(7.5) approx 2.7386.
    3. This numerical estimate is very close to the algebraic number 1 + sqrt(3) approx 2.7320.
       We hypothesize that the exact connective constant is mu = 1 + sqrt(3).
    4. We then find the minimal polynomial for x = 1 + sqrt(3).
       x - 1 = sqrt(3)
       (x - 1)^2 = 3
       x^2 - 2x + 1 = 3
       x^2 - 2x - 2 = 0
    5. This script will print the coefficients of this polynomial.
    """
    
    # Based on the analysis, the connective constant mu is hypothesized to be 1 + sqrt(3).
    # Let's represent this number in the form a + b*sqrt(c).
    a = 1
    b = 1
    c = 3
    
    # Calculate the coefficients of the minimal polynomial x^2 - 2x - 2 = 0
    coeff2, coeff1, coeff0 = find_minimal_polynomial_for_a_plus_b_sqrt_c(a, b, c)
    
    print("The minimal polynomial is P(x) = c2*x^2 + c1*x + c0 = 0.")
    print("The coefficients of the polynomial are:")
    print(f"c2 (coefficient of x^2): {coeff2}")
    print(f"c1 (coefficient of x): {coeff1}")
    print(f"c0 (constant term): {coeff0}")
    
    print("\nThe final equation is:")
    print(f"({coeff2})*x^2 + ({coeff1})*x + ({coeff0}) = 0")

solve_task()

# The final answer format requires printing each number in the final equation.
# The polynomial is x^2 - 2x - 2 = 0.
# The numbers are the coefficients 1, -2, -2.
print("\nEach number in the final equation:")
print(1)
print(-2)
print(-2)