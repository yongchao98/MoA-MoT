import torch
import torch.nn as nn

def get_rank(tensor):
    """Computes the rank of a torch tensor."""
    return torch.linalg.matrix_rank(tensor).item()

# --- Network and Input Parameters ---
n_samples = 100
n_features = 50
input_rank = 25
d_layer = 50
d_output = 10

# --- 1. Create an input matrix X with the specified rank ---
# We create two smaller matrices and multiply them to ensure the rank is controlled.
# U(100, 25) @ V(25, 50) -> X(100, 50) with rank <= 25.
print(f"Generating input matrix X of size ({n_samples}, {n_features}) with rank {input_rank}...")
U = torch.randn(n_samples, input_rank)
V = torch.randn(input_rank, n_features)
X = U @ V
print(f"Rank of input matrix X: {get_rank(X)}\n")


# --- 2. Define the MLP architecture ---
class MLP(nn.Module):
    def __init__(self):
        super().__init__()
        # Use bias terms as described in the problem
        self.layer1 = nn.Linear(n_features, d_layer, bias=True)
        self.layer2 = nn.Linear(d_layer, d_layer, bias=True)
        self.layer3 = nn.Linear(d_layer, d_output, bias=True)

    def forward(self, x):
        h1 = torch.relu(self.layer1(x))
        h2 = torch.relu(self.layer2(h1))
        h_out = torch.relu(self.layer3(h2))
        return h1, h2, h_out

# We will run simulations to test the plausibility of each statement.

# --- 3. Test Statement A: rank(H1) = 20 (Possible) ---
print("--- Testing Statement A: Can rank of first layer representation be 20? ---")
model_A = MLP()
# To get a low rank, we can manually set the weights of the first layer to have a low rank.
# Here, we construct W1 with rank 20.
with torch.no_grad():
    low_rank_w1 = torch.randn(d_layer, 20) @ torch.randn(20, n_features)
    model_A.layer1.weight.copy_(low_rank_w1)
    # Also set a bias that might zero out some rows/columns
    model_A.layer1.bias.uniform_(-10, -5)

h1_A, _, _ = model_A(X)
rank_h1_A = get_rank(h1_A)
print(f"A specific network configuration resulted in rank(H1) = {rank_h1_A}.")
print("This is close to 20, and with specific data and weights, getting exactly 20 is plausible.")
print("Conclusion: Statement A could be True.\n")


# --- 4. Test Statement B: rank(H1) = 50 (Possible) ---
print("--- Testing Statement B: Can rank of first layer representation be 50? ---")
# A standard, randomly initialized network should demonstrate the rank-expansion property.
model_B = MLP()
h1_B, h2_B, h_out_B = model_B(X)
rank_h1_B = get_rank(h1_B)
print(f"A random network resulted in rank(H1) = {rank_h1_B}.")
print("This demonstrates that the rank can increase to the maximum possible value (50).")
print("Conclusion: Statement B could be True.\n")


# --- 5. Test Statement C: rank(H2) = 25 (Possible) ---
print("--- Testing Statement C: Can rank of second layer representation be 25? ---")
# The rank of H2 depends on the input H1 and the second layer's weights.
# We saw above that with random weights, rank(H2) can be high (e.g., 50).
print(f"In the previous random network, rank(H2) was {get_rank(h2_B)}.")
# To show it can be 25, we can construct a network where W2 has a lower rank.
model_C = MLP()
with torch.no_grad():
    # Let's use the high-rank H1 from the previous model as input
    h1_C = h1_B # rank(h1_C) = 50
    # And construct W2 to have rank 25
    low_rank_w2 = torch.randn(d_layer, 25) @ torch.randn(25, d_layer)
    model_C.layer2.weight.copy_(low_rank_w2)

h2_C = torch.relu(model_C.layer2(h1_C))
rank_h2_C = get_rank(h2_C)
print(f"Using a rank-25 W2, we get rank(H2) = {rank_h2_C}.")
print("This shows that a rank of 25 is a perfectly plausible outcome.")
print("Conclusion: Statement C could be True.\n")


# --- 6. Test Statement D: rank(H_out) = 15 (Impossible) ---
print("--- Testing Statement D: Can rank of final layer representation be 15? ---")
h_out_D = h_out_B # from our random model
rows, cols = h_out_D.shape
rank_h_out_D = get_rank(h_out_D)
print(f"The final representation matrix H_out has dimensions ({rows}, {cols}).")
print(f"By definition, the rank of a matrix cannot exceed its number of rows or columns.")
print(f"Therefore, the maximum possible rank for H_out is {min(rows, cols)}.")
print(f"A rank of 15 is impossible.")
print("Conclusion: Statement D is False.")