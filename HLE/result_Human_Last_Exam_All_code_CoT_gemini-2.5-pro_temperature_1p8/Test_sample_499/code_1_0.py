import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# A helper function to get the L2 norm of a model's parameters
def get_weights_norm(model):
    norm = 0.0
    for param in model.parameters():
        norm += torch.norm(param).item() ** 2
    return np.sqrt(norm)

# A helper function to calculate the L2 distance between the parameters of two models
def get_weights_distance(model1, model2):
    dist = 0.0
    params1 = list(model1.parameters())
    params2 = list(model2.parameters())
    for i in range(len(params1)):
        dist += torch.norm(params1[i] - params2[i]).item() ** 2
    return np.sqrt(dist)
    
# 1. Define a simple model
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.layer1 = nn.Linear(1, 64)
        self.layer2 = nn.Linear(64, 64)
        self.layer3 = nn.Linear(64, 1)

    def forward(self, x):
        x = torch.relu(self.layer1(x))
        x = torch.relu(self.layer2(x))
        return self.layer3(x)

# 2. Create the data
X_train = torch.linspace(-np.pi, np.pi, 100).unsqueeze(1)
Y_train = torch.sin(X_train)

# 3. Instantiate models and set initialization scales
torch.manual_seed(42)
net_small_init = SimpleNet()
# Keep the default (Kaiming) initialization for the small model

torch.manual_seed(42)
net_large_init = SimpleNet()
# Multiply weights of the second model by a large factor
large_init_factor = 100.0
with torch.no_grad():
    for param in net_large_init.parameters():
        param.data *= large_init_factor

# Store initial states to calculate distance later
from copy import deepcopy
initial_small_net = deepcopy(net_small_init)
initial_large_net = deepcopy(net_large_init)

# 4. Training Loop
def train(model, title):
    print(f"--- Training {title} ---")
    optimizer = optim.SGD(model.parameters(), lr=0.001)
    loss_fn = nn.MSELoss()

    initial_norm = get_weights_norm(model)
    print(f"Initial L2 Norm of Weights: {initial_norm:.4f}")
    
    for epoch in range(200):
        optimizer.zero_grad()
        outputs = model(X_train)
        loss = loss_fn(outputs, Y_train)
        loss.backward()
        optimizer.step()

    final_norm = get_weights_norm(model)
    print(f"Final L2 Norm of Weights:   {final_norm:.4f}")
    return model

# Train both models
final_small_net = train(net_small_init, "Small-Init Net")
final_large_net = train(net_large_init, "Large-Init Net")
print("-" * 30)

# 5. Compare the change in weights
dist_small = get_weights_distance(initial_small_net, final_small_net)
dist_large = get_weights_distance(initial_large_net, final_large_net)

# Calculate relative distance to be fair
relative_dist_small = dist_small / get_weights_norm(initial_small_net)
relative_dist_large = dist_large / get_weights_norm(initial_large_net)

print("\n--- Results ---")
print("This demonstrates how initialization magnitude affects training dynamics.")
print("The 'lazy training' regime (approximated by perturbation theory) occurs when parameters")
print("change very little relative to their initial magnitude.\n")

print(f"Absolute change in weights (Small-Init Net): {dist_small:.4f}")
print(f"Absolute change in weights (Large-Init Net): {dist_large:.4f}")
print("\n")
print(f"Relative change in weights (Small-Init Net): {relative_dist_small:.4f}")
print(f"Relative change in weights (Large-Init Net): {relative_dist_large:.4f}")

print("\nObservation: The network with a large initialization magnitude moved its weights")
print("far less in a *relative* sense. Its behavior is dominated by its initial state,")
print("aligning with the assumptions of perturbation theory where the solution is a small")
print("correction to the initial state.")
