import math

def calculate_lower_bound():
    """
    This function calculates the constant lower bound for d(t,x).
    The lower bound is the minimum value of the function m1(u) for u in [0,1].
    The minimum is found to be at u=1.
    The function m1(u) is: (u*(3-5u) - u*sqrt(17u^2 - 22u + 9)) / 4
    """

    # The minimum of m1(u) occurs at u=1.
    u = 1.0

    # Coefficients of the polynomial inside the square root
    a = 17.0
    b = -22.0
    c = 9.0

    # Calculate each term of the equation for m1(u) at u=1
    term1 = u * (3.0 - 5.0 * u)
    discriminant = a * u**2 + b * u + c
    sqrt_term = math.sqrt(discriminant)
    numerator = term1 - u * sqrt_term
    denominator = 4.0

    # Calculate the final result
    result = numerator / denominator

    # Print the step-by-step calculation
    print("The constant lower bound is calculated by finding the minimum of m1(u) for u in [0,1].")
    print("The minimum is attained at u = 1.")
    print("\n--- Calculation Steps ---")
    print(f"The equation for the lower bound at u={u} is:")
    print(f"m1({u}) = ({u}*(3-{5}*{u}) - {u}*sqrt({a}*{u}^2 - {abs(b)}*{u} + {c})) / {denominator}")
    print(f"\n1. First term in numerator: {u} * (3 - 5*{u}) = {term1}")
    print(f"2. Term inside the square root: {a}*({u}^2) - {abs(b)}*{u} + {c} = {discriminant}")
    print(f"3. Square root term: sqrt({discriminant}) = {sqrt_term}")
    print(f"4. Numerator: {term1} - {u} * {sqrt_term} = {numerator}")
    print(f"5. Denominator: {denominator}")
    print("\n--- Final Result ---")
    print(f"Lower Bound = Numerator / Denominator = {numerator} / {denominator} = {result}")

calculate_lower_bound()
<<< -1.0 >>>