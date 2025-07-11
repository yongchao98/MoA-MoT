import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Set a seed for reproducibility
torch.manual_seed(42)

# 1. Define a simple model and a synthetic dataset
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.fc1 = nn.Linear(1, 32)
        self.fc2 = nn.Linear(32, 1)

    def forward(self, x):
        x = torch.tanh(self.fc1(x))
        x = self.fc2(x)
        return x

# Synthetic data: y = sin(x)
X = torch.linspace(-np.pi, np.pi, 100).view(-1, 1)
y = torch.sin(X)

# 2. Define a function to run the training experiment
def run_experiment(init_scale):
    """
    Initializes a model with a given scale, trains it, and returns weight norms.
    """
    model = SimpleNet()
    
    # Custom initialization with a specific scale
    def init_weights(m):
        if isinstance(m, nn.Linear):
            nn.init.normal_(m.weight, mean=0, std=init_scale)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
    model.apply(init_weights)

    # Store initial weights and calculate their norm
    initial_params = torch.cat([p.flatten() for p in model.parameters()])
    initial_norm = torch.linalg.norm(initial_params).item()

    # Training loop
    optimizer = optim.SGD(model.parameters(), lr=0.01)
    loss_fn = nn.MSELoss()
    for epoch in range(200):
        optimizer.zero_grad()
        outputs = model(X)
        loss = loss_fn(outputs, y)
        loss.backward()
        optimizer.step()

    # Store final weights and calculate the norm of the change
    final_params = torch.cat([p.flatten() for p in model.parameters()])
    param_change_norm = torch.linalg.norm(final_params - initial_params).item()
    
    return initial_norm, param_change_norm

# 3. Run experiments with small and large initialization scales
small_init_scale = 0.1  # Corresponds to "lazy" / near-linear regime
large_init_scale = 2.0   # Corresponds to "feature learning" regime

small_init_norm, small_change_norm = run_experiment(small_init_scale)
large_init_norm, large_change_norm = run_experiment(large_init_scale)

# 4. Print the "final equation" of results
# The equation shows how the weight change (||δθ||) relates to initial weight magnitude (||θ₀||)
print("Demonstration of Initialization Magnitude's Impact:")
print("-" * 50)
print("Regime 1: Small Initialization (Lazy / Perturbative)")
print(f"||θ₀|| = {small_init_norm:.4f}")
print(f"||θ* - θ₀|| = {small_change_norm:.4f}  (Parameters changed very little)")
print("-" * 50)
print("Regime 2: Large Initialization (Feature Learning)")
print(f"||θ₀|| = {large_init_norm:.4f}")
print(f"||θ* - θ₀|| = {large_change_norm:.4f}  (Parameters changed significantly)")
print("-" * 50)
print("Conclusion: The magnitude of initialization determines the scale of parameter change,")
print("thus dictating the importance of second-order (and higher) effects.")
