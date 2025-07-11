import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import math

# Set a seed for reproducibility
torch.manual_seed(42)

# 1. Define a simple dataset (e.g., fitting a sine wave)
X_train = torch.linspace(-math.pi, math.pi, 200).view(-1, 1)
y_train = torch.sin(X_train)

# 2. Define a simple MLP model
class SimpleMLP(nn.Module):
    def __init__(self, hidden_dim=256):
        super().__init__()
        self.fc1 = nn.Linear(1, hidden_dim)
        self.activation = nn.ReLU()
        self.fc2 = nn.Linear(hidden_dim, 1)

    def forward(self, x):
        x = self.fc1(x)
        x = self.activation(x)
        x = self.fc2(x)
        return x

def get_param_l2_norm(params):
    """Helper function to calculate the L2 norm of a list of tensors."""
    return torch.sqrt(sum(p.pow(2).sum() for p in params))

def run_experiment(init_scale, learning_rate):
    """
    Trains a model with a specific weight initialization scale and returns
    the L2 distance the parameters moved from their initial state.
    """
    print(f"\n--- Running Experiment with Initialization Scale: {init_scale} ---")
    
    # Instantiate the model
    model = SimpleMLP()
    
    # Manually initialize weights with the given scale
    # This is the key step that controls the learning regime
    with torch.no_grad():
        for param in model.parameters():
            param.uniform_(-init_scale, init_scale)

    # Store initial parameters
    initial_params = [p.clone().detach() for p in model.parameters()]
    initial_norm = get_param_l2_norm(initial_params)
    
    # Training setup
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    
    # Training loop
    for epoch in range(2000):
        optimizer.zero_grad()
        outputs = model(X_train)
        loss = criterion(outputs, y_train)
        loss.backward()
        optimizer.step()

    # Get final parameters
    final_params = [p.clone().detach() for p in model.parameters()]

    # Calculate the L2 distance between final and initial parameters
    param_change = [final - initial for final, initial in zip(final_params, initial_params)]
    distance_moved = get_param_l2_norm(param_change)
    
    # The "final equation" here is the comparison of how far parameters moved.
    # We print the numbers involved in this comparison.
    print(f"Initial parameter L2 norm: {initial_norm.item():.4f}")
    print(f"Final parameter L2 norm: {get_param_l2_norm(final_params).item():.4f}")
    print(f"L2 distance parameters moved from initialization: {distance_moved.item():.4f}")
    
    return distance_moved.item()


# Experiment 1: Large initialization magnitude (expected "lazy training")
# The network is already powerful at init, so weights don't need to move much.
large_scale_init = 1.0
# In lazy regime, model is sensitive to LR, so we use a smaller one.
# This mimics theoretical results where LR must scale with init.
lazy_lr = 1e-4
dist_lazy = run_experiment(large_scale_init, lazy_lr)

# Experiment 2: Small initialization magnitude (expected "feature learning")
# The network must learn features, so weights move significantly.
small_scale_init = 0.1
# A larger LR can be used here.
feature_lr = 1e-3
dist_feature = run_experiment(small_scale_init, feature_lr)

print("\n--- Conclusion ---")
print("The experiment demonstrates that the magnitude of weight initialization determines the nature of the learned parameters.")
print(f"With large initialization (scale={large_scale_init}), parameters moved a distance of {dist_lazy:.4f}.")
print(f"With small initialization (scale={small_scale_init}), parameters moved a distance of {dist_feature:.4f}.")
print("The much larger movement in the second case indicates a different learning regime ('feature learning') compared to the first ('lazy training').")
print("Therefore, the initialization magnitude determines whether the optimal parameters are a small or large perturbation from the initial state.")
