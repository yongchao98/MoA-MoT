import math

def calculate_lift_ratio(c):
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.

    The function uses the mirror image method and flow tangency conditions.
    It assumes the aerofoils are modeled as point vortices.

    Args:
        c (float): The chord length of the aerofoils.

    Returns:
        float: The calculated lift ratio L1/L2.
    """
    # Given parameters
    s = 0.5 * c  # Separation between aerofoils
    h = 0.5 * c  # Ride height

    # The problem can be solved by setting up the induced velocity equations
    # and solving for the ratio of circulations Gamma1/Gamma2, which equals L1/L2.
    #
    # The equations for induced vertical velocity at aerofoil 1 (w_ind_1) and 2 (w_ind_2) are:
    # w_ind_1 = Gamma1/(4*pi*h) + Gamma2/(2*pi) * [s/(s^2 + 4*h^2) - 1/s]
    # w_ind_2 = Gamma2/(4*pi*h) + Gamma1/(2*pi) * [1/s - s/(s^2 + 4*h^2)]
    #
    # Assuming flow tangency (w_ind_1 = w_ind_2), we can solve for Gamma1/Gamma2.
    # Let K = (1/(2*pi)) * [1/s - s/(s^2 + 4*h^2)]
    # The equation simplifies to:
    # Gamma1 * (1/(4*pi*h) - K) = Gamma2 * (1/(4*pi*h) + K)
    # Gamma1 / Gamma2 = (1/(4*pi*h) + K) / (1/(4*pi*h) - K)

    # Let's calculate the terms.
    # We can factor out 'c' as it will cancel out. Let c=1 for simplicity.
    s_norm = 0.5
    h_norm = 0.5
    
    pi = math.pi

    # Calculate the interaction term K
    K_numerator = 1/s_norm - s_norm / (s_norm**2 + 4*h_norm**2)
    K = (1 / (2 * pi)) * K_numerator

    # Calculate the self-induced ground effect term
    self_induced_term = 1 / (4 * pi * h_norm)

    # Calculate the numerator and denominator of the ratio
    ratio_numerator = self_induced_term + K
    ratio_denominator = self_induced_term - K
    
    # Calculate the final lift ratio
    lift_ratio = ratio_numerator / ratio_denominator

    # Print the steps of the final calculation
    print("This script calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.")
    print(f"Given chord c, the separation is s = 0.5c and the ride height is h = 0.5c.")
    print("The lift ratio L1/L2 is given by the formula: (A + K) / (A - K)")
    print("where A is the self-induced ground effect term and K is the interaction term.")
    print(f"\n1. Calculate the self-induced term A = 1 / (4 * pi * h):")
    print(f"   A = 1 / (4 * {pi:.4f} * {h_norm}) = {self_induced_term:.4f}")
    print(f"\n2. Calculate the interaction term K = (1 / (2 * pi)) * [1/s - s/(s^2 + 4h^2)]:")
    K_bracket_term = 1/s_norm - s_norm / (s_norm**2 + 4*h_norm**2)
    print(f"   First, the bracketed part: [1/{s_norm} - {s_norm}/({s_norm}^2 + 4*{h_norm}^2)] = {K_bracket_term:.4f}")
    print(f"   K = (1 / (2 * {pi:.4f})) * {K_bracket_term:.4f} = {K:.4f}")
    print(f"\n3. Calculate the lift ratio L1/L2 = (A + K) / (A - K):")
    print(f"   L1/L2 = ({self_induced_term:.4f} + {K:.4f}) / ({self_induced_term:.4f} - {K:.4f})")
    print(f"   L1/L2 = {ratio_numerator:.4f} / {ratio_denominator:.4f}")
    print(f"   L1/L2 = {lift_ratio:.4f}")
    
    return lift_ratio

# We can use any value for c, as it cancels out in the final ratio.
# Let's use c = 1.0 for simplicity.
final_ratio = calculate_lift_ratio(1.0)
print(f"\nFinal Answer: {final_ratio}")
