from decimal import Decimal, getcontext

def solve():
    """
    This function calculates the radius R based on the solvability condition
    of the given boundary-value problem.
    """
    # Set a high precision for calculations involving very large numbers.
    getcontext().prec = 100

    # From T = ln(10^34), we get e^T.
    e_T = Decimal('1e34')

    # The condition for the existence of a solution connects the initial values
    # (x_0, y_0, z_0) in a spherical relationship x_0^2 + y_0^2 + z_0^2 = R^2.
    # The equation for R^2 can be simplified to R^2 = 0.5 * (e^T + 1) * e^T.
    
    # Calculate R^2 using the simplified formula.
    half = Decimal('0.5')
    one = Decimal('1')
    
    R_squared = half * (e_T + one) * e_T
    
    # R is the square root of R^2.
    R = R_squared.sqrt()

    # The final equation describes the sphere of valid initial values.
    # Per the instructions, we output the numbers in this final equation.
    print("The set of initial values (x_0, y_0, z_0) must satisfy the equation:")
    # The format specifier '.5e' will present the large number in scientific notation.
    print(f"x_0^2 + y_0^2 + z_0^2 = {R_squared:.5e}")
    print("\nThis equation represents a sphere in the initial value space.")
    print(f"The radius of the sphere, R, is the square root of the value above.")
    # We print R with a high degree of precision.
    print(f"R = {R}")

solve()