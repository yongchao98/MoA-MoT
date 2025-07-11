import torch
import torch.nn as nn
import math

def demonstrate_initialization_effect(scale):
    """
    Demonstrates the effect of weight initialization scale on the initial gradient norm.

    Args:
        scale (float): The standard deviation of the normal distribution used for weight initialization.
    """
    # for reproducibility
    torch.manual_seed(42)

    # Define a simple network, a simple loss, and some random data
    n_input, n_hidden, n_output = 10, 32, 5
    batch_size = 16
    
    model = nn.Sequential(
        nn.Linear(n_input, n_hidden),
        nn.Tanh(),
        nn.Linear(n_hidden, n_output)
    )
    
    # Custom initialization
    def init_weights(m):
        if isinstance(m, nn.Linear):
            # Initialize weights from a normal distribution with the given scale (std dev)
            nn.init.normal_(m.weight, mean=0, std=scale)
            # Initialize biases to zero
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
                
    model.apply(init_weights)
    
    # Create some dummy data and target
    inputs = torch.randn(batch_size, n_input)
    targets = torch.randn(batch_size, n_output)
    
    # Calculate loss
    criterion = nn.MSELoss()
    outputs = model(inputs)
    loss = criterion(outputs, targets)
    
    # Calculate gradients
    loss.backward()
    
    # Calculate and print the L2 norm of the gradients for all parameters
    total_norm = 0
    for p in model.parameters():
        if p.grad is not None:
            param_norm = p.grad.detach().data.norm(2)
            total_norm += param_norm.item() ** 2
    total_norm = total_norm ** 0.5
    
    print(f"Initialization Scale (std dev): {scale:<10.4f} | Initial Loss: {loss.item():<10.4f} | Initial Gradient L2 Norm: {total_norm:<10.4f}")

if __name__ == '__main__':
    # He initialization for Tanh suggests std = 1/sqrt(fan_in)
    fan_in = 10
    good_scale = 1 / math.sqrt(fan_in)
    
    print("--- Demonstrating the effect of weight initialization scale on initial gradients ---")
    # A very small scale (likely to cause vanishing gradients)
    demonstrate_initialization_effect(scale=0.01)
    
    # A "good" scale based on common heuristics (like Xavier/He)
    demonstrate_initialization_effect(scale=good_scale)
    
    # A very large scale (likely to cause exploding gradients)
    demonstrate_initialization_effect(scale=1.0)
    print("---------------------------------------------------------------------------------")
    print("Notice how the gradient norm changes dramatically with the initialization scale.")
    print("This initial gradient, along with the Hessian (curvature), is what a perturbation theory analysis would study.")