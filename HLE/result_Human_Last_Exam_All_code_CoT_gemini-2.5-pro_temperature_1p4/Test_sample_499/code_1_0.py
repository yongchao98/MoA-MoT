import torch
import torch.nn as nn

# Define a simple feedforward neural network
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        self.fc1 = nn.Linear(10, 50)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(50, 2)

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x

# --- Case 1: Standard (Small) Initialization ---
net_small_init = SimpleNet()
# PyTorch default initialization (Kaiming uniform) is relatively small.

# --- Case 2: Large Initialization ---
net_large_init = SimpleNet()
def init_weights_large(m):
    if isinstance(m, nn.Linear):
        # Initialize weights with a large standard deviation
        nn.init.normal_(m.weight, mean=0, std=10.0)
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)
net_large_init.apply(init_weights_large)


# Create some dummy data
input_tensor = torch.randn(1, 10)
target = torch.tensor([1], dtype=torch.long)
criterion = nn.CrossEntropyLoss()

# --- Calculate gradients for the small initialization network ---
output_small = net_small_init(input_tensor)
loss_small = criterion(output_small, target)
loss_small.backward()

total_grad_norm_small = 0
for param in net_small_init.parameters():
    if param.grad is not None:
        total_grad_norm_small += param.grad.data.norm(2).item() ** 2
total_grad_norm_small = total_grad_norm_small ** 0.5


# --- Calculate gradients for the large initialization network ---
output_large = net_large_init(input_tensor)
loss_large = criterion(output_large, target)
loss_large.backward()

total_grad_norm_large = 0
for param in net_large_init.parameters():
    if param.grad is not None:
        total_grad_norm_large += param.grad.data.norm(2).item() ** 2
total_grad_norm_large = total_grad_norm_large ** 0.5

print(f"Total gradient norm with small (standard) initialization: {total_grad_norm_small:.4f}")
print(f"Total gradient norm with large initialization: {total_grad_norm_large:.4f}")
print("\nThis demonstrates that the magnitude of weight initialization directly impacts the gradient,")
print("which is the first-order term in the perturbation theory expansion of the loss function.")
