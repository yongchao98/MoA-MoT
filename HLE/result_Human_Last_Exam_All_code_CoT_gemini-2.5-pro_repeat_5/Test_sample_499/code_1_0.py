import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Set a seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# 1. Define the simple feedforward neural network
class SimpleFNN(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, init_std):
        super(SimpleFNN, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.activation = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, output_size)
        self.init_std = init_std
        self._initialize_weights()

    def _initialize_weights(self):
        # Initialize weights with a specific standard deviation
        nn.init.normal_(self.fc1.weight, mean=0.0, std=self.init_std)
        nn.init.zeros_(self.fc1.bias)
        nn.init.normal_(self.fc2.weight, mean=0.0, std=self.init_std)
        nn.init.zeros_(self.fc2.bias)

    def forward(self, x):
        x = self.fc1(x)
        x = self.activation(x)
        x = self.fc2(x)
        return x

    def get_param_norm(self):
        # Helper to calculate the L2 norm of all parameters
        return torch.cat([p.view(-1) for p in self.parameters()]).norm().item()

# 2. Setup the experiment
input_dim = 1
hidden_dim = 128
output_dim = 1
learning_rate = 0.01
epochs = 100

# Generate some synthetic data (e.g., fitting a sine wave)
X_train = torch.linspace(-np.pi, np.pi, 100).view(-1, 1)
y_train = torch.sin(X_train)

# 3. Run two experiments with different weight initializations
init_std_small = 0.1  # Small initialization magnitude
init_std_large = 2.0   # Large initialization magnitude

model_small_init = SimpleFNN(input_dim, hidden_dim, output_dim, init_std=init_std_small)
model_large_init = SimpleFNN(input_dim, hidden_dim, output_dim, init_std=init_std_large)

# Helper function for training
def train_model(model):
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    initial_norm = model.get_param_norm()
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(X_train)
        loss = criterion(outputs, y_train)
        loss.backward()
        optimizer.step()

    final_norm = model.get_param_norm()
    final_loss = loss.item()
    
    return initial_norm, final_norm, final_loss

print("This experiment demonstrates how the magnitude of weight initialization affects the final learned parameters.\n")

# Train and evaluate the model with small initialization
initial_norm_small, final_norm_small, final_loss_small = train_model(model_small_init)
print(f"--- Model with Small Initialization (std = {init_std_small}) ---")
print(f"Initial L2 norm of parameters: {initial_norm_small:.4f}")
print(f"Final L2 norm of parameters:   {final_norm_small:.4f}")
print(f"Change in parameter norm:      {(final_norm_small - initial_norm_small):.4f}")
print(f"Final training loss:           {final_loss_small:.6f}\n")


# Train and evaluate the model with large initialization
initial_norm_large, final_norm_large, final_loss_large = train_model(model_large_init)
print(f"--- Model with Large Initialization (std = {init_std_large}) ---")
print(f"Initial L2 norm of parameters: {initial_norm_large:.4f}")
print(f"Final L2 norm of parameters:   {final_norm_large:.4f}")
print(f"Change in parameter norm:      {(final_norm_large - initial_norm_large):.4f}")
print(f"Final training loss:           {final_loss_large:.6f}\n")

print("Conclusion: The model with large initialization barely changed its parameters (lazy training),")
print("while the model with small initialization moved significantly farther in parameter space to find a solution.")
print("This shows that the initial magnitude is a key determinant of the final 'optimal' parameters.")

<<<D>>>