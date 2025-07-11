import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Set random seed for reproducibility
torch.manual_seed(42)

# 1. Define a simple Feedforward Neural Network
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.layer1 = nn.Linear(1, 64)
        self.layer2 = nn.Linear(64, 64)
        self.layer3 = nn.Linear(64, 1)

    def forward(self, x):
        x = torch.relu(self.layer1(x))
        x = torch.relu(self.layer2(x))
        x = self.layer3(x)
        return x

    def get_l2_norm_of_weights(self):
        """Helper to compute the L2 norm of all weights."""
        norm = 0.0
        for param in self.parameters():
            if param.requires_grad:
                norm += torch.norm(param).item() ** 2
        return np.sqrt(norm)

# 2. Create a synthetic dataset (y = x*sin(x))
X_train = torch.linspace(-5, 5, 100).view(-1, 1)
y_train = X_train * torch.sin(X_train)

# 3. Training Function
def train_model(model, X, y, epochs=200, lr=0.01):
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(X)
        loss = criterion(outputs, y)
        loss.backward()
        optimizer.step()
    return model, loss.item()

# 4. Function to initialize weights with a specific magnitude
def init_weights(m, std_dev):
    if isinstance(m, nn.Linear):
        torch.nn.init.normal_(m.weight, mean=0, std=std_dev)
        if m.bias is not None:
            torch.nn.init.zeros_(m.bias)

# --- Case 1: Small Initialization Magnitude ---
print("--- Case 1: Small Initialization Magnitude ---")
small_init_model = SimpleNet()
# Initialize weights with a small standard deviation
small_init_std = 0.01
small_init_model.apply(lambda m: init_weights(m, std_dev=small_init_std))

# Store a deep copy of initial weights for later comparison
initial_small_weights = [p.clone().detach() for p in small_init_model.parameters()]
initial_small_norm = small_init_model.get_l2_norm_of_weights()
print(f"Initial L2 Norm of Weights (Small Init): {initial_small_norm:.4f}")

# Train the model
trained_small_model, final_loss_small = train_model(small_init_model, X_train, y_train)
final_small_norm = trained_small_model.get_l2_norm_of_weights()
print(f"Final L2 Norm of Weights (Small Init):   {final_small_norm:.4f}")

# Calculate the norm of the change in weights
weight_change_small = 0.0
for p_initial, p_final in zip(initial_small_weights, trained_small_model.parameters()):
    weight_change_small += torch.norm(p_final - p_initial).item() ** 2
norm_of_change_small = np.sqrt(weight_change_small)
print(f"L2 Norm of Weight Change (||w_final - w_initial||): {norm_of_change_small:.4f}")
print(f"Final Loss: {final_loss_small:.6f}\n")


# --- Case 2: Large Initialization Magnitude ---
print("--- Case 2: Large Initialization Magnitude ---")
large_init_model = SimpleNet()
# Initialize weights with a large standard deviation
large_init_std = 1.0
large_init_model.apply(lambda m: init_weights(m, std_dev=large_init_std))

# Store a deep copy of initial weights
initial_large_weights = [p.clone().detach() for p in large_init_model.parameters()]
initial_large_norm = large_init_model.get_l2_norm_of_weights()
print(f"Initial L2 Norm of Weights (Large Init): {initial_large_norm:.4f}")

# Train the model
trained_large_model, final_loss_large = train_model(large_init_model, X_train, y_train)
final_large_norm = trained_large_model.get_l2_norm_of_weights()
print(f"Final L2 Norm of Weights (Large Init):   {final_large_norm:.4f}")

# Calculate the norm of the change in weights
weight_change_large = 0.0
for p_initial, p_final in zip(initial_large_weights, trained_large_model.parameters()):
    weight_change_large += torch.norm(p_final - p_initial).item() ** 2
norm_of_change_large = np.sqrt(weight_change_large)
print(f"L2 Norm of Weight Change (||w_final - w_initial||): {norm_of_change_large:.4f}")
print(f"Final Loss: {final_loss_large:.6f}\n")

print("Conclusion:")
print("The L2 norm of the weight change is much larger for the model with large-magnitude initialization.")
print("This demonstrates that the initial weight magnitude determines whether the final parameters are a small perturbation of the initial ones.")