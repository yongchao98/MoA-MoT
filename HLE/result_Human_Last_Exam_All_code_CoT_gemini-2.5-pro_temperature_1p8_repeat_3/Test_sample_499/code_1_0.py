import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# 1. Define a simple Neural Network
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.layer1 = nn.Linear(1, 32)
        self.activation = nn.Tanh()
        self.layer2 = nn.Linear(32, 1)

    def forward(self, x):
        x = self.layer1(x)
        x = self.activation(x)
        x = self.layer2(x)
        return x

# Function to initialize weights with a specific magnitude (std dev)
def init_weights(m, magnitude):
    if isinstance(m, nn.Linear):
        nn.init.normal_(m.weight, mean=0, std=magnitude)
        nn.init.constant_(m.bias, 0)

# Function to calculate the L2 norm of model parameters
def get_param_norm(model):
    return torch.cat([p.view(-1) for p in model.parameters()]).norm().item()

# 2. Setup the training problem
X_train = torch.linspace(-np.pi, np.pi, 100).view(-1, 1)
y_train = torch.sin(X_train)
loss_fn = nn.MSELoss()

# 3. Experiment with different initialization magnitudes
magnitudes = {
    "Small Magnitude (σ=0.1)": 0.1,
    "Large Magnitude (σ=1.0)": 1.0,
}

print("Demonstrating how initialization magnitude affects the final learned parameters.\n")

for name, mag in magnitudes.items():
    print(f"--- Running experiment with: {name} ---")

    # Create a new model instance for each experiment
    model = SimpleNet()
    
    # Apply the specific weight initialization
    model.apply(lambda m: init_weights(m, magnitude=mag))
    
    initial_norm = get_param_norm(model)
    print(f"Initial L2 norm of weights: {initial_norm:.4f}")
    
    # Train the model
    optimizer = optim.SGD(model.parameters(), lr=0.01)
    for epoch in range(500):
        optimizer.zero_grad()
        y_pred = model(X_train)
        loss = loss_fn(y_pred, y_train)
        loss.backward()
        optimizer.step()

    final_norm = get_param_norm(model)
    print(f"Final L2 norm of weights after training: {final_norm:.4f}")
    print(f"Final loss: {loss.item():.6f}\n")

# 4. Final Conclusion based on the theory
print("="*60)
print("Conclusion:")
print("The experiments show that starting from different weight magnitudes,")
print("even with everything else being identical (architecture, data, optimizer),")
print("results in different final sets of parameters (as shown by the different final norms).")
print("\nIn a perturbation theory interpretation of neural networks, the analysis often")
print("expands the network function in a series around zero weights. The magnitude")
print("of the initial weights acts as the perturbation scale parameter.")
print("\nThis scale parameter determines the relative influence of first-order effects")
print("(akin to kernel methods) versus second-order and higher effects (feature learning).")
print("Therefore, the property that determines the optimal parameters is:")
print("\nD. magnitude of weight initialization")
print("="*60)