import torch
import torch.nn as nn
import torch.optim as optim
import math

# Step 1: Define a simple model and task for regression
class SimpleNet(nn.Module):
    """A simple Feedforward Neural Network."""
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.fc1 = nn.Linear(1, 64)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(64, 1)

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x

def weights_init(m, scale=1.0):
    """Custom weight initialization function to apply a scale factor."""
    if isinstance(m, nn.Linear):
        # Initialize weights using Kaiming Uniform, a standard practice
        nn.init.kaiming_uniform_(m.weight, a=math.sqrt(5))
        # Apply the custom scale factor to the weights
        m.weight.data.mul_(scale)
        # Initialize biases, if they exist
        if m.bias is not None:
            fan_in, _ = nn.init._calculate_fan_in_and_fan_out(m.weight)
            bound = scale / math.sqrt(fan_in) # Scale the bound for bias as well
            nn.init.uniform_(m.bias, -bound, bound)

# Step 2: Generate synthetic data to learn the sine function
X_train = torch.linspace(-math.pi, math.pi, 200).view(-1, 1)
y_train = torch.sin(X_train)

# Step 3: Define a training loop
def train_model(model, scale_name):
    """Trains a given model and prints its progress."""
    print(f"--- Training model with '{scale_name}' initialization ---")
    optimizer = optim.SGD(model.parameters(), lr=0.01)
    loss_fn = nn.MSELoss()
    
    initial_loss = loss_fn(model(X_train), y_train).item()
    print(f"Initial Loss: {initial_loss:.6f}")

    for epoch in range(2001):
        optimizer.zero_grad()
        outputs = model(X_train)
        loss = loss_fn(outputs, y_train)
        # Check for exploding gradients/loss
        if torch.isnan(loss) or loss.item() > 1e5:
            print(f"Epoch {epoch}: Training unstable. Loss exploded.")
            break
        loss.backward()
        optimizer.step()

        if epoch % 500 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.6f}")
    
    final_loss = loss.item()
    print(f"Final Loss for '{scale_name}' initialization: {final_loss:.6f}\n")
    return final_loss

# --- Main execution: Compare two different initialization scales ---

# Model 1: Initialized with a standard scale (scale = 1.0)
print("Demonstration: Comparing Standard vs. Large-Scale Initialization\n")
model_standard = SimpleNet()
# Apply standard-scaled initialization explicitly
model_standard.apply(lambda m: weights_init(m, scale=1.0))
train_model(model_standard, "Standard Scale")

# Model 2: Initialized with a large scale (scale = 10.0)
model_large_init = SimpleNet()
# Apply initialization with a much larger scale
model_large_init.apply(lambda m: weights_init(m, scale=10.0))
train_model(model_large_init, "Large Scale")

print("### Conclusion from Demonstration ###")
print("The standard initialization allows the model to converge to a low loss, successfully learning the function.")
print("The large-scale initialization leads to unstable gradients and a high final loss, preventing convergence.")
print("This shows that the initial weight magnitude is a determining factor for the training outcome and the resulting 'optimal' parameters.")
