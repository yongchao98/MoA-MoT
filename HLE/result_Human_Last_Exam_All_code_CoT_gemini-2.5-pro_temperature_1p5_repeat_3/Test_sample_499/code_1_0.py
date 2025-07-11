import numpy as np

def run_experiment():
    """
    Demonstrates the effect of initialization magnitude on final parameters.
    """
    # A simple "neural network": f(x) = w2 * max(0, w1 * x)
    # We want to learn a target function, y = 2.0 * x
    # An ideal solution is w1=1.0, w2=2.0 (or any w1=c, w2=2.0/c for c>0)

    # Training data
    X_train = np.array([[1.0], [2.0], [3.0]])
    Y_train = np.array([[2.0], [4.0], [6.0]])

    # Hyperparameters
    learning_rate = 0.01
    epochs = 100
    
    # We fix the random seed to ensure the 'direction' of the initial
    # random vector is the same, only its magnitude changes.
    np.random.seed(0)
    w1_base = np.random.randn()
    w2_base = np.random.randn()

    def train(w1_init, w2_init, case_name):
        """Train a simple model and print the results."""
        w1, w2 = w1_init, w2_init
        
        initial_norm = np.sqrt(w1**2 + w2**2)
        print(f"--- {case_name} ---")
        print(f"Initial Weights: w1 = {w1:.4f}, w2 = {w2:.4f}")
        print(f"Initial Parameter Norm: {initial_norm:.4f}")

        for epoch in range(epochs):
            # Forward pass
            h = X_train * w1
            h_relu = np.maximum(h, 0)
            y_pred = h_relu * w2
            
            # Loss
            loss = np.mean((y_pred - Y_train)**2)

            # Backward pass (Gradient calculation)
            grad_y_pred = 2.0 * (y_pred - Y_train) / len(Y_train)
            grad_w2 = np.sum(grad_y_pred * h_relu)
            grad_h_relu = grad_y_pred * w2
            grad_h = grad_h_relu.copy()
            grad_h[h < 0] = 0
            grad_w1 = np.sum(grad_h * X_train)

            # Update weights
            w1 -= learning_rate * grad_w1
            w2 -= learning_rate * grad_w2
        
        final_norm = np.sqrt(w1**2 + w2**2)
        param_change = np.sqrt((w1 - w1_init)**2 + (w2 - w2_init)**2)
        
        print(f"Final Weights:   w1 = {w1:.4f}, w2 = {w2:.4f}")
        print(f"Final Loss: {loss:.6f}")
        print(f"Distance moved in parameter space: {param_change:.4f}")
        print(f"Ratio of distance moved to initial norm: {param_change / initial_norm:.4f}\n")

    # Case 1: Small Initialization Magnitude (feature learning regime)
    train(w1_base * 0.1, w2_base * 0.1, "Small Initialization Magnitude")

    # Case 2: Large Initialization Magnitude (lazy/kernel regime)
    train(w1_base * 10.0, w2_base * 10.0, "Large Initialization Magnitude")

# Run the demonstration
run_experiment()