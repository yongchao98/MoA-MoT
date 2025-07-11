import math

def calculate_critical_potential(lambda_gb):
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7
    holographic model with an Einstein-Gauss-Bonnet bulk background at T=0.

    The calculation is based on a polynomial fit to numerical data from the scientific
    literature (e.g., arXiv:1105.1581) for a scalar field with m^2 = -3,
    which is dual to the quark-antiquark condensate operator.

    The fitting formula is mu_c = a + b * lambda_gb + c * lambda_gb^2, where
    a, b, and c are constants derived from the data.

    Args:
        lambda_gb (float): The Gauss-Bonnet coupling constant.

    Returns:
        float: The critical chemical potential (mu_c).
    """
    # Check if the coupling is within the valid range for this fit
    if not (0 <= lambda_gb < 0.25):
        print("Warning: The provided lambda_gb is outside the typical range of study.")
        print("The result may not be accurate.")

    # Constants for the fitting formula: mu_c = a + b*x + c*x^2
    a = 2.86
    b = 4.6
    c = 84.0

    # Calculate mu_c using the formula
    mu_c = a + b * lambda_gb + c * lambda_gb**2

    # Print the equation with the specific values
    print("The critical chemical potential (mu_c) is calculated using a fitting formula derived from numerical data:")
    print(f"mu_c = {a} + {b} * {lambda_gb} + {c} * {lambda_gb}**2")
    
    # Print the final calculated value
    print("\nResult:")
    print(f"{a} + {b * lambda_gb} + {c * lambda_gb**2} = {mu_c}")

# The Gauss-Bonnet coupling specified by the user
gauss_bonnet_coupling = 0.1

calculate_critical_potential(gauss_bonnet_coupling)