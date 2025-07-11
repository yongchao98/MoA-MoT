import numpy as np
from scipy.special import lambertw

def solve_droplet_charge():
    """
    Calculates the total charge on a spherical droplet based on the provided equations.
    """
    # Step 1: Define the given constants
    sigma_0 = 7.43e-7  # Units: e/nm
    R_0 = 30.0         # Units: nm
    q_i = 2 * np.pi    # Constant

    # Step 2: Calculate the necessary values involving the Lambert W function.
    # The problem defines omega = W(1), the Lambert W function evaluated at 1.
    # The lambertw function returns a complex number, so we take the real part.
    omega = lambertw(1).real

    # The simplified integral requires calculating Y_pi = W(exp(q_i * pi)).
    Y_pi = lambertw(np.exp(q_i * np.pi)).real

    # Step 3: Use the final simplified formula for Q.
    # The derivation leads to the expression:
    # Q = (2 * pi * sigma_0 * R_0 / q_i^2) * [ (Y_pi * (ln(Y_pi) - 1)) / (1 + Y_pi) + omega ]
    
    # Step 4: Print the final equation with all symbolic and numerical values.
    print("Based on the simplification, the final expression for the total charge Q is:")
    print("Q = (2 * pi * sigma_0 * R_0 / q_i^2) * [ (Y_pi * (ln(Y_pi) - 1)) / (1 + Y_pi) + omega ]")
    print("\nWith the following values:")
    print(f"sigma_0 = {sigma_0} e/nm")
    print(f"R_0 = {R_0} nm")
    print(f"q_i = {q_i:.6f}")
    print(f"omega = W(1) = {omega:.6f}")
    print(f"Y_pi = W(exp(q_i * pi)) = {Y_pi:.6f}")

    print("\nSubstituting the numerical values into the equation gives:")
    ln_Y_pi = np.log(Y_pi)
    
    # We explicitly show each number in the final equation string as requested.
    # To maintain clarity, we show the results of the sub-expressions.
    term1_num_val = Y_pi * (ln_Y_pi - 1)
    term1_den_val = 1 + Y_pi
    term1_val = term1_num_val / term1_den_val
    bracket_val = term1_val + omega
    prefactor_val = (2 * np.pi * sigma_0 * R_0) / (q_i**2)

    print(f"Q = (2 * {np.pi:.6f} * {sigma_0} * {R_0} / {q_i:.6f}**2) * [ ({Y_pi:.6f} * (ln({Y_pi:.6f}) - 1)) / (1 + {Y_pi:.6f}) + {omega:.6f} ]")
    print(f"Q = ({prefactor_val:.6e}) * [ (({term1_num_val:.6f}) / ({term1_den_val:.6f})) + {omega:.6f} ]")
    print(f"Q = ({prefactor_val:.6e}) * [ {term1_val:.6f} + {omega:.6f} ]")
    print(f"Q = ({prefactor_val:.6e}) * [ {bracket_val:.6f} ]")

    # Step 5: Calculate the final numerical result for Q.
    Q = prefactor_val * bracket_val

    print(f"\nThe calculated total charge on the droplet is:")
    print(f"Q = {Q:.3e} e")
    
    return Q

# Execute the function and store the result
total_charge = solve_droplet_charge()
print(f"\n<<<{total_charge}>>>")
