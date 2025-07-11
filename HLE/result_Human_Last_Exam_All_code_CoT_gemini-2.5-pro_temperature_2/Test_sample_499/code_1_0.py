import math

def analyze_network_stability():
    """
    Analyzes a key relationship from the perturbation theory of neural networks.

    In this theory, a critical parameter, let's call it `chi`, determines
    if signals propagate stably through the network. `chi > 1` leads to a
    "chaotic" regime with exploding gradients and a very rough loss landscape.
    `chi < 1` leads to an "ordered" regime with vanishing gradients.
    The optimal state for training is at the "edge of chaos", where `chi = 1`.

    A simplified form of the equation for `chi` is:
    chi = C_activation * sigma_w^2

    where:
    - C_activation is a constant derived from the activation function (e.g., for tanh, C is approx 1).
    - sigma_w is the standard deviation of the weight initialization distribution.

    This shows that `sigma_w`, which represents the MAGNITUDE of the weight initialization,
    is the parameter that directly determines this critical property.
    """

    # For many common activation functions like tanh, the contribution is ~1.
    C_activation = 1.0

    # The critical value for sigma_w is then the one that makes chi = 1.
    # We choose sigma_w = 1.0 to be at the critical point.
    sigma_w = 1.0

    # Calculate the chaos parameter 'chi'
    chi = C_activation * (sigma_w ** 2)

    # Print the explanation and the final equation with all numbers.
    print("This script illustrates how weight magnitude determines network stability.")
    print("The goal is to be at the 'edge of chaos', where the parameter chi = 1.")
    print("The governing equation is: chi = C_activation * (sigma_w ** 2)")
    print("\nLet's plug in the numbers for the critical state:")
    # We use f-string formatting to print the equation with its values.
    # Note that we are printing each individual number that goes into the calculation.
    print(f"{chi:.1f} = {C_activation:.1f} * ({sigma_w:.1f} ** 2)")
    print("\nThis demonstrates that the weight initialization magnitude (represented by sigma_w) directly controls the network's dynamical regime, which is a key finding of the second-order perturbation analysis.")

# Execute the analysis.
analyze_network_stability()