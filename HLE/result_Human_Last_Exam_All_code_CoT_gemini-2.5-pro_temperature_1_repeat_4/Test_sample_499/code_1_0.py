import sympy

def demonstrate_weight_magnitude_effect():
    """
    This function demonstrates how the magnitude of weight initialization
    determines the second-order properties (curvature) of the loss landscape.
    """
    # 1. Define symbolic variables for a simple 2-layer linear network
    # y_pred = w2 * w1 * x
    w1, w2, x, y_true = sympy.symbols('w1 w2 x y_true')

    # Define the network's prediction
    y_pred = w2 * w1 * x

    # 2. Define a simple loss function (Mean Squared Error)
    loss = (y_true - y_pred)**2
    print("--- Symbolic Analysis ---")
    print(f"Simple Network Prediction: y_pred = {y_pred}")
    print(f"Loss Function: L = {loss}\n")

    # 3. Calculate the second derivatives of the loss with respect to the weights.
    # These form the Hessian matrix, which represents the landscape's curvature.
    H11 = sympy.diff(loss, w1, 2)
    H22 = sympy.diff(loss, w2, 2)
    H12 = sympy.diff(loss, w1, w2)

    print("Second derivative of Loss w.r.t. w1 (Hessian component H11):")
    print(f"d²L/dw1² = {H11}")
    print("\nSecond derivative of Loss w.r.t. w2 (Hessian component H22):")
    print(f"d²L/dw2² = {H22}\n")
    print("These equations show that the curvature of the loss landscape (the values in the")
    print("Hessian matrix) is directly proportional to the square of the other weights' magnitudes.\n")

    # 4. Substitute numerical values to show the effect of initialization magnitude.
    print("--- Numerical Demonstration ---")
    print("Let's evaluate the curvature for a sample data point (x=2, y_true=5).\n")
    # Define a substitution dictionary for the data point
    subs_data = {x: 2, y_true: 5}

    # Case 1: Small weight initialization magnitude
    w1_small, w2_small = 0.1, 0.1
    subs_small_weights = {w1: w1_small, w2: w2_small}
    H11_val_small = H11.subs(subs_data).subs(subs_small_weights)
    H22_val_small = H22.subs(subs_data).subs(subs_small_weights)
    print("Case 1: Small Weight Initialization (w1=0.1, w2=0.1)")
    print(f"Curvature w.r.t. w1 (H11) = 2 * ({w2_small})**2 * ({subs_data[x]})**2 = {H11_val_small}")
    print(f"Curvature w.r.t. w2 (H22) = 2 * ({w1_small})**2 * ({subs_data[x]})**2 = {H22_val_small}")
    print("Result: The curvature is very small, indicating a flat loss landscape.\n")


    # Case 2: Large weight initialization magnitude
    w1_large, w2_large = 10.0, 10.0
    subs_large_weights = {w1: w1_large, w2: w2_large}
    H11_val_large = H11.subs(subs_data).subs(subs_large_weights)
    H22_val_large = H22.subs(subs_data).subs(subs_large_weights)
    print("Case 2: Large Weight Initialization (w1=10.0, w2=10.0)")
    print(f"Curvature w.r.t. w1 (H11) = 2 * ({w2_large})**2 * ({subs_data[x]})**2 = {H11_val_large}")
    print(f"Curvature w.r.t. w2 (H22) = 2 * ({w1_large})**2 * ({subs_data[x]})**2 = {H22_val_large}")
    print("Result: The curvature is very large, indicating a sharp and chaotic loss landscape.\n")

    print("Conclusion: The magnitude of weight initialization fundamentally determines the second-order")
    print("properties of the network, which is the central idea in a perturbation theory interpretation.")


if __name__ == '__main__':
    demonstrate_weight_magnitude_effect()

<<<D>>>