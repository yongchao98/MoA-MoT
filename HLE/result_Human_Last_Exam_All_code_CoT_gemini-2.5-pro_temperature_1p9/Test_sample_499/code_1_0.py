import torch
import torch.nn as nn
import numpy as np

# Define a simple feedforward neural network
class SimpleNet(nn.Module):
    def __init__(self, hidden_size=64):
        super(SimpleNet, self).__init__()
        self.layer1 = nn.Linear(1, hidden_size)
        self.activation = nn.ReLU()
        self.layer2 = nn.Linear(hidden_size, 1)

    def forward(self, x):
        x = self.layer1(x)
        x = self.activation(x)
        x = self.layer2(x)
        return x

# Custom weight initialization function based on a specific standard deviation (magnitude)
def init_weights_with_magnitude(module, std_dev):
    if isinstance(module, nn.Linear):
        # Initialize weights from a normal distribution with the specified magnitude
        torch.nn.init.normal_(module.weight, mean=0.0, std=std_dev)
        # Initialize biases to zero for consistency
        if module.bias is not None:
            torch.nn.init.constant_(module.bias, 0)

# --- Simulation Parameters ---
# We will test different magnitudes (standard deviations) for weight initialization
init_magnitudes = [0.1, 1.0, 5.0]
# Other parameters are fixed
learning_rate = 0.01
num_epochs = 200
# Fix random seeds for reproducible results across runs
torch.manual_seed(42)
np.random.seed(42)

# --- Create Synthetic Data (fitting a sine wave) ---
X_train = torch.linspace(-np.pi, np.pi, 200).view(-1, 1)
Y_train = torch.sin(X_train)

print("Investigating the effect of weight initialization magnitude on final model performance.")
print(f"Fixed Parameters: Learning Rate = {learning_rate}, Epochs = {num_epochs}\n")

# --- Run experiment for each initialization magnitude ---
for magnitude in init_magnitudes:
    # 1. Create a new model instance for each experiment to start fresh
    model = SimpleNet()
    
    # 2. Apply the custom weight initialization with the specified magnitude
    model.apply(lambda m: init_weights_with_magnitude(m, std_dev=magnitude))
    
    # 3. Define the loss function and optimizer
    criterion = nn.MSELoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
    
    # 4. Train the model
    for epoch in range(num_epochs):
        outputs = model(X_train)
        loss = criterion(outputs, Y_train)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    final_loss = loss.item()
    
    # 5. Report the result for the current initialization magnitude
    print(f"Initialization Magnitude (std dev): {magnitude:<4.1f} -> Final Loss (MSE): {final_loss:.6f}")

print("\nConclusion: The final performance (and thus the learned 'optimal' parameters) is")
print("highly dependent on the initial weight magnitude, even when all other factors are constant.")
print("This aligns with perturbation theory, where the initial state determines the landscape to be explored.")

<<<D>>>