import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet based on the derived formula.
    """
    # Given constants
    sigma_0 = 7.43e-7  # in units of e/nm
    R_0 = 30.0         # in units of nm

    # Calculate omega = W(1), the Lambert W function evaluated at 1
    # lambertw returns a complex number, so we take the real part.
    omega = np.real(lambertw(1))

    # Calculate the constant C from the simplified expression
    C = omega / (1 + omega)**3

    # Calculate pi^3
    pi_cubed = np.pi**3

    # Calculate the total charge Q using the derived formula: Q = C * pi^3 * sigma_0 * R_0
    Q = C * pi_cubed * sigma_0 * R_0

    # Print the breakdown of the calculation
    print("Calculating the total charge Q using the formula:")
    print("Q = (omega / (1 + omega)**3) * pi**3 * sigma_0 * R_0")
    print("\nSubstituting values:")
    print(f"omega = W(1) = {omega}")
    print(f"C = omega / (1 + omega)**3 = {C}")
    print(f"pi^3 = {pi_cubed}")
    print(f"sigma_0 = {sigma_0} e/nm")
    print(f"R_0 = {R_0} nm")
    print(f"\nFinal calculation:")
    print(f"Q = {C} * {pi_cubed} * {sigma_0} * {R_0}")
    
    # Print the final result
    print(f"\nTotal charge Q = {Q} e")

if __name__ == "__main__":
    calculate_total_charge()
    # The final numerical result required for the answer box
    # Q = 0.0001018903543265778
    # Let's provide this in the requested format
    # The prompt doesn't specify precision, so using a reasonable number of sig figs
    # For example, 1.019e-4
    # The user didn't ask for a specific format, I will just output the number.
    # The instructions say "directly return the answer with the format <<<answer content>>> at the end of your response"
    # So I will calculate the value and put it in that format.
    sigma_0 = 7.43e-7
    R_0 = 30.0
    omega = np.real(lambertw(1))
    C = omega / (1 + omega)**3
    pi_cubed = np.pi**3
    Q = C * pi_cubed * sigma_0 * R_0
    print(f"\n<<<{Q}>>>")
