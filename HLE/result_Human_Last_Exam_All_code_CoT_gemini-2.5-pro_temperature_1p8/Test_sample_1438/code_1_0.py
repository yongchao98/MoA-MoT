def berezin_integral(c0, c1):
    """
    Simulates the Berezin integral for a function of a single Grassmann variable.
    The integrand is f(eta) = c0 + c1*eta.
    The integral is defined as integral(f(eta) deta) = c1.

    Args:
        c0 (float or int): The constant term in the function.
        c1 (float or int): The coefficient of the Grassmann variable 'eta'.

    Returns:
        float or int: The result of the integration.
    """
    # The Berezin integral picks out the coefficient of the linear term.
    return c1

def main():
    """
    Main function to demonstrate the value of the Grassmann integral measure.
    """
    print("In the path integral formalism for fermions, integration over a Grassmann variable 'eta' is defined by rules.")
    print("An arbitrary function f(eta) can be written as f(eta) = c0 + c1*eta because eta^2 = 0.")
    print("The Berezin integral is defined as: integral( (c0 + c1*eta) deta ) = c1.\n")

    # Case 1: Integral of a constant, f(eta) = 1 = 1 + 0*eta
    c0_const = 1
    c1_const = 0
    result_const = berezin_integral(c0_const, c1_const)
    print(f"Let's calculate the integral of a constant, f(eta) = {c0_const}:")
    print(f"integral( {c0_const} deta ) = {result_const}")

    # Case 2: The integral that normalizes the measure, f(eta) = eta = 0 + 1*eta
    c0_var = 0
    c1_var = 1
    result_var = berezin_integral(c0_var, c1_var)
    print(f"\nNow, let's calculate the integral that defines the measure's normalization, f(eta) = eta:")
    # The final equation prints each component number
    print(f"integral( ({c0_var} + {c1_var}*eta) deta ) = {result_var}")
    print("\nThis result establishes the normalization of the measure. The value which maintains the Pauli principle is 1 in this context.")

if __name__ == "__main__":
    main()
