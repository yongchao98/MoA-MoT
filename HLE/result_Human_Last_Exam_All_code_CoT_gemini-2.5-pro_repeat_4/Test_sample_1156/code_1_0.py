import sympy

def find_normalized_density():
    """
    This function calculates the normalized density for the given map.
    It starts with the known unnormalized density, calculates the normalization
    constant by integrating it over the interval [0, 1], and then
    constructs the final expression for the normalized density.
    """
    # Define the symbolic variable
    x = sympy.symbols('x')

    # The unnormalized invariant density from literature
    a_un = 2
    b_un = -1
    rho_unnormalized = a_un / sympy.sqrt(x) + b_un

    # Calculate the normalization constant C
    try:
        C = sympy.integrate(rho_unnormalized, (x, 0, 1))
    except (sympy.SympifyError, TypeError, ValueError) as e:
        print(f"Error during integration: {e}")
        return

    # Calculate the coefficients of the normalized density
    a_norm = sympy.Rational(a_un, C)
    b_norm = sympy.Rational(b_un, C)

    # Output the final equation for the normalized density
    print("The normalised density of the invariant measure is rho(x) = a/sqrt(x) + b")
    print("where the coefficients are:")
    print(f"a = {a_un}/{C} = {a_norm}")
    print(f"b = {b_un}/{C} = {b_norm}")
    print("\nThe final equation is:")
    print(f"rho(x) = ({a_norm})/sqrt(x) + ({b_norm})")

find_normalized_density()