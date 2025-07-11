import math

def find_minimal_polynomial_from_squared_length():
    """
    Calculates the minimal polynomial for a length L, given L^2 = a + b*sqrt(5).

    For this problem, the shortest geodesic distance L for an ant starting and ending
    at the same vertex of a dodecahedron (without passing through other vertices)
    has a squared length L^2 = 7 - 3*sqrt(5).
    """
    # L^2 = a + b*sqrt(5)
    a = 7
    b = -3

    # We are looking for the minimal polynomial of L.
    # Let x = L, so x^2 = a + b*sqrt(5).
    # To find the polynomial with rational coefficients, we isolate the sqrt term and square it.
    # (x^2 - a)^2 = (b*sqrt(5))^2
    # x^4 - 2*a*x^2 + a^2 = b^2 * 5
    # x^4 - 2*a*x^2 + (a^2 - 5*b^2) = 0

    # Calculate the coefficients of the polynomial: c4*x^4 + c2*x^2 + c0 = 0
    c4 = 1
    c2 = -2 * a
    c0 = a**2 - 5 * b**2

    # The minimal polynomial is x^4 - 14x^2 + 4.
    # We will print the final equation and its coefficients.
    
    print(f"The squared length of the shortest path is L^2 = {a} + {b}*sqrt(5).")
    print("The minimal polynomial p(x) for the length L satisfies p(L) = 0.")
    
    # Build the string for the equation
    equation_str = f"{c4}*x^4 + ({c2})*x^2 + ({c0}) = 0"
    print(f"\nThe final equation is: {equation_str}")

    print("\nThe numbers in the final equation are:")
    print(c4)
    print(c2)
    print(c0)

# Run the function
find_minimal_polynomial_from_squared_length()