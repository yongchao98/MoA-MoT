import math

def demonstrate_perturbation_theory_role():
    """
    This function conceptually demonstrates why weight initialization magnitude
    is key in the perturbation theory view of neural networks.
    """
    
    # In perturbation theory (e.g., "lazy training" or NTK regime),
    # we assume the change in weights (delta_w) is small compared to
    # the initial weights (w_0).
    
    # Let's model a weight update with gradient descent: delta_w = -learning_rate * gradient
    # The crucial factor is the *relative* change: delta_w / w_0
    
    learning_rate = 0.01
    # Let's assume a representative gradient value for this thought experiment
    gradient = 5.0
    
    # Calculate the absolute change in weight, which is constant in this example
    delta_w = -learning_rate * gradient
    
    print("--- Perturbation Theory & Initialization Magnitude ---")
    print(f"Goal: Analyze the network's behavior by treating training as a small perturbation.")
    print(f"Condition: The relative weight change |Δw / w₀| must be small.")
    print(f"Fixed learning rate = {learning_rate}, Fixed gradient = {gradient}")
    print(f"The absolute weight change (Δw) is {delta_w}")
    print("-" * 50)
    
    # --- Case 1: Large Initialization Magnitude (triggers "lazy" / perturbative regime) ---
    initial_weight_magnitude_large = 10.0
    
    relative_change_large_init = abs(delta_w / initial_weight_magnitude_large)
    
    print("Case 1: Large Weight Initialization (Perturbative Regime)")
    print("Conceptual Equation: Relative_Change = |Δw / w₀|")
    print("Plugging in the numbers:")
    # Outputting each number in the final equation as requested
    print(f"Relative_Change = |{delta_w} / {initial_weight_magnitude_large}| = {relative_change_large_init:.4f}")
    print("Result: The relative change is small. The perturbation theory view is valid.")
    print("The network's optimal parameters are found on a linearized loss landscape.")

    print("-" * 50)
    
    # --- Case 2: Small Initialization Magnitude (triggers "rich" / feature learning regime) ---
    initial_weight_magnitude_small = 0.1
    
    relative_change_small_init = abs(delta_w / initial_weight_magnitude_small)
    
    print("Case 2: Small Weight Initialization (Feature Learning Regime)")
    print("Conceptual Equation: Relative_Change = |Δw / w₀|")
    print("Plugging in the numbers:")
    # Outputting each number in the final equation as requested
    print(f"Relative_Change = |{delta_w} / {initial_weight_magnitude_small}| = {relative_change_small_init:.4f}")
    print("Result: The relative change is large. The simple perturbation theory view breaks down.")
    print("The network learns complex features, changing its internal function significantly.")

    print("\nConclusion: The magnitude of weight initialization determines which regime the network is in,")
    print("and thus the nature of its optimal parameters under a perturbative analysis.")

demonstrate_perturbation_theory_role()
<<<D>>>