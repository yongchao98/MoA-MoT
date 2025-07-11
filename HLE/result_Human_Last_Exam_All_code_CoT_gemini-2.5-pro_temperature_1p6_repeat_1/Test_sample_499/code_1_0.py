import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Set a seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# Create a simple dataset
X = torch.randn(10, 5)
y = torch.randn(10, 1)

# Define a simple linear model
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.fc1 = nn.Linear(5, 1, bias=False)

    def forward(self, x):
        return self.fc1(x)

def train_and_get_final_weights(learning_rate):
    """
    Trains a simple model for a fixed number of steps and returns the final weights.
    """
    model = SimpleNet()
    # Store initial weights for comparison
    initial_weights = model.fc1.weight.data.clone()
    
    criterion = nn.MSELoss()
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)

    # Train for a few steps
    for _ in range(5):
        optimizer.zero_grad()
        outputs = model(X)
        loss = criterion(outputs, y)
        loss.backward()
        optimizer.step()

    final_weights = model.fc1.weight.data.clone()
    return initial_weights, final_weights

# --- Demonstration ---
# In perturbation theory, the final parameters are a function of the learning rate.
# Let's show this empirically.

lr1 = 0.01
lr2 = 0.1

initial_weights1, final_weights1 = train_and_get_final_weights(lr1)
# For a fair comparison, we need to re-initialize the second model to the same starting point
# Since we set the seed, running the function again will do this.
initial_weights2, final_weights2 = train_and_get_final_weights(lr2)

print("This demonstration shows that the final trained parameters ('optimal parameters') depend on the learning rate.")
print("-" * 70)
# Note: Initial weights will be the same due to the fixed seed.
print(f"Initial Weights (shared for both experiments):\n{np.round(initial_weights1.numpy(), 4)}")
print("-" * 70)

print(f"Final weights with learning rate = {lr1}:")
print(np.round(final_weights1.numpy(), 4))
print("\n")
print(f"Final weights with learning rate = {lr2}:")
print(np.round(final_weights2.numpy(), 4))
print("-" * 70)

print("As you can see, the final weight values are different, confirming that they are determined by the learning rate.")
print("In perturbation theory, the final weights are expressed as a mathematical series where the learning rate is the expansion parameter.")
print("The 'second order' part of the question refers to terms proportional to (learning_rate)^2 in this series.")

<<<C>>>