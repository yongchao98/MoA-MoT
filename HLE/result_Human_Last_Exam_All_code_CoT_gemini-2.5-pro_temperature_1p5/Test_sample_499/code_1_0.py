import torch
import torch.nn as nn
import numpy as np
import warnings

# Suppress user warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

def create_mlp(input_dim, hidden_dim, output_dim, depth):
    """Creates a simple Multi-Layer Perceptron (MLP)."""
    layers = [nn.Linear(input_dim, hidden_dim), nn.ReLU()]
    for _ in range(depth - 1):
        layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.ReLU()])
    layers.append(nn.Linear(hidden_dim, output_dim))
    return nn.Sequential(*layers)

def analyze_network_hessian(model, model_name, input_dim):
    """Computes and analyzes the Hessian of the network's loss."""
    print(f"--- Analyzing {model_name} ---")
    
    # We need to flatten the model's parameters into a single vector for the Hessian calculation
    flat_params, unflatten_fn = torch.func.flatten_params(model.parameters())
    num_params = len(flat_params)
    print(f"Total number of parameters: {num_params}")

    # Dummy data for loss calculation
    x_sample = torch.randn(1, input_dim)
    y_sample = torch.randn(1, 1)

    # A stateless loss function that takes a flat parameter vector
    def flattened_loss_func(p_vec):
        # Unflatten the vector back into the model's parameter structure
        unflattened_params = unflatten_fn(p_vec)
        # Use a stateless version of the model to compute output
        # `functional_call` replaces the model's parameters with the provided ones for a single forward pass
        y_pred = torch.func.functional_call(model, {name: p for name, p in zip(model.state_dict(), unflattened_params)}, (x_sample,))
        return nn.MSELoss()(y_pred, y_sample)

    # Compute the 2D Hessian matrix of the loss function w.r.t. the flattened parameters
    hessian_2d = torch.autograd.functional.hessian(flattened_loss_func, flat_params)
    
    # Calculate eigenvalues of the Hessian
    try:
        eigenvalues = np.linalg.eigvalsh(hessian_2d.detach().numpy())
        # We look at the magnitude of eigenvalues
        eigenvalues_abs = np.abs(eigenvalues)
        
        # To avoid division by zero for singular matrices, we consider a small epsilon
        min_eig_mag = np.min(eigenvalues_abs[eigenvalues_abs > 1e-9])
        max_eig_mag = np.max(eigenvalues_abs)
        
        if min_eig_mag > 1e-9:
            condition_number = max_eig_mag / min_eig_mag
            print(f"Max eigenvalue (magnitude): {max_eig_mag:.4f}")
            print(f"Min eigenvalue (magnitude): {min_eig_mag:.4f}")
            print(f"Hessian Condition Number: {condition_number:.2e}\n")
        else:
            print("Hessian is singular (min non-zero eigenvalue is too small). Condition number is effectively infinite.\n")
            
    except np.linalg.LinAlgError:
        print("Could not compute eigenvalues for the Hessian.\n")

# --- Main analysis ---
torch.manual_seed(42)
INPUT_DIM = 10

# Case 1: Wide, Shallow Network (small depth-to-width ratio)
depth_wide = 3
width_wide = 100
ratio_wide = depth_wide / width_wide
print(f"Creating WIDE network (Depth={depth_wide}, Width={width_wide}, Ratio={ratio_wide:.3f})")
wide_model = create_mlp(input_dim=INPUT_DIM, hidden_dim=width_wide, output_dim=1, depth=depth_wide)
analyze_network_hessian(wide_model, "Wide Network", INPUT_DIM)

# Case 2: Deep, Narrow Network (large depth-to-width ratio)
depth_deep = 10
width_deep = 15
ratio_deep = depth_deep / width_deep
print(f"Creating DEEP network (Depth={depth_deep}, Width={width_deep}, Ratio={ratio_deep:.3f})")
deep_model = create_mlp(input_dim=INPUT_DIM, hidden_dim=width_deep, output_dim=1, depth=depth_deep)
analyze_network_hessian(deep_model, "Deep Network", INPUT_DIM)