import torch

def analyze_initialization_effect():
    """
    Demonstrates how the magnitude of weight initialization affects the
    predicted optimal parameters under a 2nd-order perturbation theory.
    """
    # Define a simple problem: given x_true, predict y_true
    x_true = torch.tensor([2.0])
    y_true = torch.tensor([0.9])
    
    # The model is y_pred = tanh(w * x + b)
    # The loss is L = (y_pred - y_true)^2

    def compute_and_print_analysis(w_init, b_init):
        """
        Computes gradient, Hessian, and predicted optimal parameters
        for a given initialization point (w_init, b_init).
        """
        # Define parameters at the initial point
        p_init = torch.tensor([w_init, b_init])
        
        # Define the loss function in terms of a parameter vector p=[w, b]
        def loss_fn(p):
            w, b = p[0], p[1]
            y_pred = torch.tanh(w * x_true + b)
            return (y_pred - y_true)**2

        # --- Perturbation Theory Calculations ---
        # Equation for optimal parameters: p* = p_init - H(p_init)⁻¹ * g(p_init)
        
        # 1. Calculate the gradient g(p_init)
        # torch.autograd.functional.jacobian computes the gradient for a scalar function
        g = torch.autograd.functional.jacobian(loss_fn, p_init)
        
        # 2. Calculate the Hessian H(p_init)
        H = torch.autograd.functional.hessian(loss_fn, p_init)

        print(f"Initial Parameters (p_init): w={p_init[0].item():.4f}, b={p_init[1].item():.4f}")
        print(f"Gradient at init (g): [w_grad={g[0].item():.4f}, b_grad={g[1].item():.4f}]")
        print(f"Hessian at init (H):")
        print(f"  [[{H[0,0].item():.4f}, {H[0,1].item():.4f}],")
        print(f"   [{H[1,0].item():.4f}, {H[1,1].item():.4f}]]")

        # 3. Calculate the optimal parameters p*
        try:
            # Solve H * dp = -g for the update dp = p* - p_init
            dp = -torch.linalg.solve(H, g)
            p_optimal = p_init + dp
            print(f"Predicted Optimal Parameters (p*): w*={p_optimal[0].item():.4f}, b*={p_optimal[1].item():.4f}")
        except torch.linalg.LinAlgError:
            print("Hessian is singular. Cannot compute optimal parameters via inversion.")


    # --- Case 1: Small Magnitude Initialization ---
    print("--- Case 1: Small Magnitude Weight Initialization ---")
    # Here, tanh is in its linear region.
    compute_and_print_analysis(w_init=0.1, b_init=0.0)

    # --- Case 2: Large Magnitude Initialization ---
    print("\n--- Case 2: Large Magnitude Weight Initialization ---")
    # Here, tanh is saturated, leading to tiny gradients (vanishing gradient problem).
    compute_and_print_analysis(w_init=3.0, b_init=0.0)
    
    print("\nAs shown, changing the initial weight magnitude dramatically alters the gradient and Hessian,")
    print("which in turn changes the predicted optimal parameters derived from the perturbation theory.")

if __name__ == '__main__':
    analyze_initialization_effect()