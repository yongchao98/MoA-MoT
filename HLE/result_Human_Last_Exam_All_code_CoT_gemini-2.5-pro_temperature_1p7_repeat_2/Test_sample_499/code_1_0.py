import torch
from torch.func import grad, hessian

def second_order_analysis(w_init, x, y):
    """
    Performs a second-order perturbation analysis to find optimal weights.

    Args:
        w_init (torch.Tensor): The initial weights of the network.
        x (torch.Tensor): The input data.
        y (torch.Tensor): The target data.
    """
    print(f"\n--- Analysis for Initial Weights: w_init = {w_init.tolist()} ---")

    # Define the loss function in terms of the weights for a simple model:
    # y_pred = w2 * tanh(w1 * x)
    def compute_loss(weights):
        w1, w2 = weights[0], weights[1]
        y_pred = w2 * torch.tanh(w1 * x)
        loss = (y_pred - y) ** 2
        return loss

    # 1. Compute the Gradient (g) and Hessian (H) at the initial weights w_init
    g = grad(compute_loss)(w_init)
    H = hessian(compute_loss)(w_init)

    # 2. Compute the inverse of the Hessian.
    # Add a small value (regularization) to the diagonal for numerical stability.
    try:
        H_inv = torch.inverse(H)
    except torch.linalg.LinAlgError:
        print("Hessian is singular. Adding regularization (1e-6 * I) for inversion.")
        H_inv = torch.inverse(H + torch.eye(H.shape[0]) * 1e-6)

    # 3. Calculate the parameter update using the formula: Δw = -H⁻¹ * g
    # Unsqueeze g to make it a column vector for matrix multiplication
    delta_w = -H_inv @ g.unsqueeze(1)
    delta_w = delta_w.squeeze(1) # Squeeze back to a 1D tensor

    # 4. Calculate the predicted optimal parameters: w_opt = w_init + Δw
    w_opt = w_init + delta_w

    # 5. Print the full equation with all its components
    print("The optimal parameters (w_opt) are estimated using the formula: w_opt = w_init - H_inv @ g")
    print("\nComponent Values:")
    print(f"  w_init (Initial Parameters)  = {w_init.tolist()}")
    print(f"  g (Gradient at w_init)     = {g.tolist()}")
    print(f"  H_inv (Inverse Hessian)      = \n{H_inv.tolist()}")
    print("\nFinal Equation:")
    print(f"  w_opt {w_opt.tolist()} = {w_init.tolist()} - [Inverse Hessian] @ {g.tolist()}")


if __name__ == '__main__':
    # Define a single data point for our simple model
    x_data = torch.tensor(1.5)
    y_data = torch.tensor(0.8)

    # Case 1: Initialize with a small magnitude
    w_init_small = torch.tensor([0.5, 0.5])
    second_order_analysis(w_init_small, x_data, y_data)

    # Case 2: Initialize with a larger magnitude
    w_init_large = torch.tensor([2.0, 2.0])
    second_order_analysis(w_init_large, x_data, y_data)

    print("\nConclusion: As shown by the different results for 'w_opt',")
    print("the magnitude of weight initialization is the property that determines")
    print("the optimal parameters under this second-order perturbation interpretation.")
