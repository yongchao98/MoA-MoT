import math

def explain_perturbation_theory_and_weights():
    """
    Explains why weight initialization magnitude is key for the perturbation theory
    view of neural networks.
    """

    print("### Perturbation Theory and Optimal Parameters in FNNs ###")
    print("\nThe perturbation theory interpretation of a neural network analyzes its training dynamics by expanding the loss function as a Taylor series around the initial weights (w₀).")
    print("Let L(w) be the loss for weights w. The expansion is:")
    print("L(w) ≈ L(w₀) + ∇L(w₀)ᵀ(w - w₀) + 0.5 * (w - w₀)ᵀ H(w₀) (w - w₀)\n")

    print("This approximation is valid when the final weights 'w' remain close to the initial weights 'w₀'. The property that most directly controls this is the magnitude (or scale) of the initial weights.\n")

    # Conceptual demonstration
    initialization_scale_alpha = 10.0  # A large scale factor
    width = 512
    # Standard initialization is often ~ 1/sqrt(fan_in)
    standard_init_magnitude = 1 / math.sqrt(width)
    large_init_magnitude = initialization_scale_alpha * standard_init_magnitude

    print(f"For a network of width {width}:")
    print(f"- A standard initialization magnitude might be: {standard_init_magnitude:.4f}")
    print(f"- A large initialization magnitude (alpha = {initialization_scale_alpha}) would be: {large_init_magnitude:.4f}\n")

    print("The regime the network operates in is determined by this magnitude:")
    print("1. Large Magnitude -> 'Lazy' or 'Perturbative' Regime: Optimal parameters are a small perturbation from initial parameters.")
    print("2. Small Magnitude -> 'Rich' or 'Feature Learning' Regime: Optimal parameters may be far from initial ones.\n")

    # The required "equation" format
    print("Therefore, we can write a conceptual equation:")
    final_weights = "w_final"
    initial_weights = "w_0"
    perturbation_term = "delta_w"
    magnitude = large_init_magnitude

    # We print each number/variable in the final equation
    print(f"Under a large initialization magnitude ({magnitude:.4f}), the final parameters can be expressed as:")
    print(f"{final_weights} = {initial_weights} + {perturbation_term}, where '{perturbation_term}' is small.")
    print("\nThis confirms that the magnitude of weight initialization (D) is the determining property.")

# Execute the explanation
explain_perturbation_theory_and_weights()