import torch

# Set a seed for reproducibility of random initializations
torch.manual_seed(42)

# --- Plan ---
# We will demonstrate that the magnitude of weight initialization determines the importance
# of the second-order term (Hessian) in a network's perturbation expansion.
#
# A neural network's output f(w) near its initial weights w0 can be approximated by:
# f(w0 + dw) ≈ f(w0) + [Gradient]⋅dw + 0.5⋅dwᵀ⋅[Hessian]⋅dw
#
# We will show that the ratio of ||Hessian|| / ||Gradient|| changes based on the
# initial weight scale. A higher ratio implies the second-order term is more
# significant for determining the network's behavior and its optimal parameters.

def analyze_network_taylor_terms(init_scale, input_dim=10, hidden_dim=20, output_dim=1):
    """
    Initializes a network with a given weight scale, and computes the norms
    of the Gradient and Hessian of its output function with respect to the weights.
    """
    # 1. Define and initialize a simple network
    # We define the function explicitly: f(x) = w2 @ tanh(w1 @ x)
    w1 = torch.randn(hidden_dim, input_dim) * init_scale
    w2 = torch.randn(output_dim, hidden_dim) * init_scale
    w1.requires_grad = True
    w2.requires_grad = True
    params = [w1, w2]
    num_params = w1.numel() + w2.numel()

    # 2. Define a sample input vector
    x = torch.randn(input_dim, 1)

    # 3. Compute the network's output
    output = (w2 @ torch.tanh(w1 @ x)).squeeze()

    # 4. Compute the Gradient (related to the first-order term)
    # This is the Jacobian of the scalar output with respect to all weight parameters.
    grad_tuple = torch.autograd.grad(output, params, create_graph=True)
    flat_grad = torch.cat([g.view(-1) for g in grad_tuple])
    grad_norm = torch.linalg.norm(flat_grad)

    # 5. Compute the Hessian (the second-order term)
    # We build the Hessian matrix by taking the gradient of each element of the gradient vector.
    hessian = torch.zeros(num_params, num_params)
    for i in range(num_params):
        # Calculate the gradient of the i-th element of the (already computed) flat_grad
        grad_of_grad_i_tuple = torch.autograd.grad(flat_grad[i], params, retain_graph=True)
        # Flatten and store it as a row in the Hessian matrix
        hessian[i] = torch.cat([g.view(-1) for g in grad_of_grad_i_tuple])
    hessian_norm = torch.linalg.norm(hessian)

    return grad_norm, hessian_norm

# --- Analysis Execution ---
print("Analyzing how weight initialization magnitude affects a network's Taylor expansion terms.")
print("-" * 80)
print("This code demonstrates that initialization scale controls the importance of the")
print("second-order (Hessian) term, which is central to a second-order perturbation theory.")
print("-" * 80)

# Case 1: Large initialization magnitude (approaching the "lazy" or linear regime)
large_scale = 5.0
grad_norm_large, hessian_norm_large = analyze_network_taylor_terms(large_scale)
ratio_large = hessian_norm_large / grad_norm_large

print(f"For a LARGE initialization scale ({large_scale}):")
print(f"  - Norm of Gradient (1st order): {grad_norm_large:.4f}")
print(f"  - Norm of Hessian (2nd order): {hessian_norm_large:.4f}")
print(f"  - Ratio (Hessian/Gradient Norm): {ratio_large:.4f}")
print("-" * 80)

# Case 2: Small initialization magnitude (the "rich" or non-linear regime)
small_scale = 0.1
grad_norm_small, hessian_norm_small = analyze_network_taylor_terms(small_scale)
ratio_small = hessian_norm_small / grad_norm_small

print(f"For a SMALL initialization scale ({small_scale}):")
print(f"  - Norm of Gradient (1st order): {grad_norm_small:.4f}")
print(f"  - Norm of Hessian (2nd order): {hessian_norm_small:.4f}")
print(f"  - Ratio (Hessian/Gradient Norm): {ratio_small:.4f}")
print("-" * 80)

print("\nConclusion:")
print(f"The ratio for the small scale ({ratio_small:.2f}) is larger than for the large scale ({ratio_large:.2f}).")
print("This shows that for smaller weight initializations, the function has more significant")
print("curvature (Hessian) relative to its slope (Gradient). Therefore, a second-order")
print("analysis is more critical, and the nature of the optimal parameters is determined by")
print("this non-linear behavior, which is governed by the magnitude of weight initialization.")
