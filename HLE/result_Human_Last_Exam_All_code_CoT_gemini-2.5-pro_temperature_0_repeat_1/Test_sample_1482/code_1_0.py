import math

def generate_it3_mf_formulation():
    """
    Generates and prints the mathematical formulation for the upper bound
    of a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF).

    The formulation uses a Gaussian paradigm where the center and standard deviation
    are linear functions of the secondary variable 'u'.
    """
    # Define parameters for the linear functions that describe the Gaussian's parameters.
    # For the center: c(u) = c1*u + c0
    c1 = 0.5
    c0 = 5.0

    # For the standard deviation of the upper bound: σ_upper(u) = sigma1*u + sigma0
    sigma1 = 0.2
    sigma0 = 1.0

    # Construct the final equation string using an f-string.
    # This equation represents the upper bound of the membership function, μ_upper(x, u).
    equation = (
        f"μ_upper(x, u) = exp( - ( (x - ({c1}*u + {c0}))^2 / "
        f"(2 * ({sigma1}*u + {sigma0})^2) ) )"
    )

    # --- Output Section ---
    print("This script formulates the equation for a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF).")
    print("\n--- CONTEXT ---")
    print("The formulation describes the upper bound of the membership grade, μ_upper(x, u),")
    print("for a fixed primary input 'x' and a secondary input 'u'.")
    print("It is based on a Gaussian function whose center 'c' and standard deviation 'σ' are")
    print("themselves linear functions of the secondary variable 'u'.")

    print("\n--- PARAMETERS ---")
    print(f"Center function: c(u) = {c1}*u + {c0}")
    print(f"Standard Deviation function: σ(u) = {sigma1}*u + {sigma0}")

    print("\n--- FINAL FORMULATION ---")
    print("The resulting mathematical formulation is:")
    print(equation)


if __name__ == "__main__":
    generate_it3_mf_formulation()
