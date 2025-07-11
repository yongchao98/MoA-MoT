import numpy as np

def demonstrate_lipschitz_effect():
    """
    Demonstrates how the Lipschitz constant of an activation function
    governs its sensitivity to perturbations.
    """
    print("--- Analyzing the role of the Lipschitz Constant in Perturbation Theory ---")
    print("Perturbation theory studies how a function's output changes when its input is slightly perturbed.")
    print("The Lipschitz constant (L) provides an upper bound on this change: |f(x_2) - f(x_1)| <= L * |x_2 - x_1|.")
    print("A lower Lipschitz constant implies more stability under perturbation, a key feature for optimal parameters.\n")

    # --- Define inputs and perturbation ---
    x = 0.5
    delta = 0.01
    x_perturbed = x + delta

    # --- Case 1: Activation with a low Lipschitz constant (L=1) ---
    print("--- Case 1: tanh Activation (Lipschitz Constant L = 1.0) ---")
    activation_tanh = np.tanh
    y = activation_tanh(x)
    y_perturbed = activation_tanh(x_perturbed)
    output_change_tanh = np.abs(y_perturbed - y)
    
    # The inequality we are checking is |f(x+delta) - f(x)| <= L * |delta|
    bound_tanh = 1.0 * np.abs(delta)

    print(f"Original input x: {x}")
    print(f"Perturbation delta: {delta}")
    print(f"Perturbed input x': {x_perturbed:.4f}")
    print("\n--- The Perturbation Equation: |f(x') - f(x)| <= L * |delta| ---")
    # Outputting each "number in the equation" as requested
    print(f"Actual Output Change |f(x') - f(x)|: {output_change_tanh:.6f}")
    print(f"Theoretical Bound L * |delta|:       {bound_tanh:.6f}")
    print(f"Is change <= bound? {output_change_tanh <= bound_tanh}")
    print("Result: With a low L, the output change is small and well-controlled.\n")


    # --- Case 2: Activation with a high Lipschitz constant (L=5) ---
    alpha = 5.0
    print(f"--- Case 2: Steep Leaky ReLU (Lipschitz Constant L = {alpha}) ---")
    def steep_leaky_relu(val):
        return val if val > 0 else alpha * val

    y = steep_leaky_relu(x)
    y_perturbed = steep_leaky_relu(x_perturbed)
    # Using a negative input to see the effect of the steep negative slope
    x_neg = -0.5
    x_neg_perturbed = x_neg + delta
    y_neg = steep_leaky_relu(x_neg)
    y_neg_perturbed = steep_leaky_relu(x_neg_perturbed)
    output_change_leaky = np.abs(y_neg_perturbed - y_neg)
    bound_leaky = alpha * np.abs(delta)

    print(f"Original input x: {x_neg}")
    print(f"Perturbation delta: {delta}")
    print(f"Perturbed input x': {x_neg_perturbed:.4f}")
    print("\n--- The Perturbation Equation: |f(x') - f(x)| <= L * |delta| ---")
    # Outputting each "number in the equation" as requested
    print(f"Actual Output Change |f(x') - f(x)|: {output_change_leaky:.6f}")
    print(f"Theoretical Bound L * |delta|:       {bound_leaky:.6f}")
    print(f"Is change <= bound? {output_change_leaky <= bound_leaky}")
    print("Result: With a high L, the same input perturbation causes a much larger output change.\n")
    
    print("--- Conclusion ---")
    print("The Lipschitz constant is a fundamental property that directly controls the network's stability against perturbations.")
    print("Therefore, it is a critical parameter in any perturbation theory analysis aimed at finding 'optimal' (i.e., stable and robust) network parameters.")


if __name__ == '__main__':
    demonstrate_lipschitz_effect()
<<<G>>>