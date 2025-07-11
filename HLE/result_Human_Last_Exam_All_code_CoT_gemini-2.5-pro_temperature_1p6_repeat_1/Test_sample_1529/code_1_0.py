# First, ensure you have the cypari2 library installed:
# pip install cypari2

from cypari2 import Pari

def solve_curve_discriminant():
    """
    Finds the minimal discriminant of the curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    pari = Pari()

    # The curve is defined by y^2 = f(x), where f(x) is the polynomial below.
    # f(x) = 1*x^6 + 0*x^5 + 0*x^4 + 2*x^3 + 4*x^2 + 4*x + 1
    # We provide the coefficients of f(x) to PARI as a vector, from the highest degree to the constant term.
    f_coeffs = [1, 0, 0, 2, 4, 4, 1]

    # ellinit finds the elliptic curve associated with the hyperelliptic curve
    E = pari.ellinit(f_coeffs)

    # ellminimalmodel finds the minimal model of the curve E
    E_minimal = pari.ellminimalmodel(E)

    # Extract the coefficients and invariants of the minimal model
    # The model is y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    a_coeffs = E_minimal.a
    a1, a2, a3, a4, a6 = a_coeffs[0], a_coeffs[1], a_coeffs[2], a_coeffs[3], a_coeffs[4]

    # The b-invariants are used to calculate the discriminant
    b2 = E_minimal.b2
    b4 = E_minimal.b4
    b6 = E_minimal.b6
    b8 = E_minimal.b8

    # The minimal discriminant
    minimal_discriminant = E_minimal.disc
    
    print(f"The minimal model of the elliptic curve is given by the equation:")
    print(f"y^2 + ({a1})*x*y + ({a3})*y = x^3 + ({a2})*x^2 + ({a4})*x + ({a6})")
    
    print("\nCalculating the minimal discriminant:")
    print(f"The b-invariants are: b2={b2}, b4={b4}, b6={b6}, b8={b8}")
    
    # Display the formula for the discriminant with the numbers plugged in
    print("The discriminant Delta is calculated as: Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6")
    print(f"Delta = -({b2})^2 * ({b8}) - 8*({b4})^3 - 27*({b6})^2 + 9*({b2})*({b4})*({b6})")
    
    # Show the calculation step by step
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    print(f"Delta = ({term1}) + ({term2}) + ({term3}) + ({term4})")
    
    print(f"\nThe minimal discriminant of the curve is: {minimal_discriminant}")

solve_curve_discriminant()