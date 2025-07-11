import numpy as np

def demonstrate_perturbation_effect():
    """
    This function demonstrates how the Lipschitz constant of an activation function
    affects a network's stability under weight perturbation.
    """
    print("--- Perturbation Theory Demonstration ---")
    print("We will compare two activation functions to see how their 'steepness' (Lipschitz constant)")
    print("affects the output when a weight is slightly perturbed.\n")

    # --- Setup ---
    # An input value
    x = 1.0
    # An initial weight
    w = 0.5
    # A small perturbation to the weight
    delta = 0.01
    w_perturbed = w + delta

    # --- Activation Functions ---
    # 1. Standard tanh function. Its Lipschitz constant is 1.
    def activation_low_lipschitz(val):
        return np.tanh(val)

    # 2. A "steeper" tanh function. Its Lipschitz constant is 3.
    def activation_high_lipschitz(val):
        return np.tanh(3 * val)

    # --- Neuron Simulation: Low Lipschitz (Smoother) ---
    print("Case 1: Neuron with low Lipschitz constant activation (tanh(x))")
    output_initial_low = activation_low_lipschitz(w * x)
    output_perturbed_low = activation_low_lipschitz(w_perturbed * x)
    change_low = np.abs(output_perturbed_low - output_initial_low)

    print(f"Initial Weight (w): {w}")
    print(f"Perturbed Weight (w + delta): {w_perturbed}")
    print(f"Initial Output: {output_initial_low:.6f}")
    print(f"Perturbed Output: {output_perturbed_low:.6f}")
    print("Final Equation: |f(w_perturbed) - f(w)| = |{:.6f} - {:.6f}| = {:.6f}\n".format(
        output_perturbed_low, output_initial_low, change_low))

    # --- Neuron Simulation: High Lipschitz (Steeper) ---
    print("Case 2: Neuron with high Lipschitz constant activation (tanh(3x))")
    output_initial_high = activation_high_lipschitz(w * x)
    output_perturbed_high = activation_high_lipschitz(w_perturbed * x)
    change_high = np.abs(output_perturbed_high - output_initial_high)

    print(f"Initial Weight (w): {w}")
    print(f"Perturbed Weight (w + delta): {w_perturbed}")
    print(f"Initial Output: {output_initial_high:.6f}")
    print(f"Perturbed Output: {output_perturbed_high:.6f}")
    print("Final Equation: |f(w_perturbed) - f(w)| = |{:.6f} - {:.6f}| = {:.6f}\n".format(
        output_perturbed_high, output_initial_high, change_high))

    print("--- Conclusion ---")
    print(f"The change in output for the smoother function ({change_low:.6f}) is much smaller")
    print(f"than for the steeper function ({change_high:.6f}).")
    print("This shows that a lower Lipschitz constant leads to greater stability under perturbation,")
    print("a key property for finding robust optimal parameters in a second-order analysis.")

if __name__ == '__main__':
    demonstrate_perturbation_effect()
<<<G>>>