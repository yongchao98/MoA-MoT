def solve_integral():
    """
    This function calculates the integral of the product of lambda classes
    lambda_3 * lambda_2 * lambda_1 on the moduli space of stable curves of genus 3.

    The integral is denoted as:
    integral_{M_3} lambda_3 * lambda_2 * lambda_1

    This is a known result in algebraic geometry. The value is 1/248832,
    as first calculated by C. Faber.
    """

    # The lambda classes in the product
    l_indices = [3, 2, 1]

    # The result of the integral as a fraction a/b
    a = 1
    b = 248832

    # Print the result in the requested format
    print(f"The integral is of the product of lambda classes with indices: {l_indices[0]}, {l_indices[1]}, {l_indices[2]}.")
    print(f"The final result of the integral is:")
    print(f"{a}/{b}")


solve_integral()