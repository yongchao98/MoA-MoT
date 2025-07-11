import numpy as np

# --- Step 1: Define activation functions and their derivatives ---
# Tanh has a Lipschitz constant of 1, as its derivative is max 1.
def tanh(x):
    return np.tanh(x)

def tanh_derivative(x):
    return 1 - np.tanh(x)**2

# A scaled Tanh has a larger Lipschitz constant.
# If f(x) = k*tanh(x), its derivative is k*(1-tanh(x)^2), so L=k.
def scaled_tanh(x, k=2.0):
    return k * np.tanh(x)

def scaled_tanh_derivative(x, k=2.0):
    return k * (1 - np.tanh(x)**2)

# --- Step 2: Define a simple 1-hidden-layer neural network ---
# f(x) = w2 * activation(w1 * x)
def simple_nn(x, w1, w2, activation_fn):
    return w2 * activation_fn(w1 * x)

# --- Step 3: Setup parameters and a small perturbation ---
x_input = 1.0
w1_initial = 0.5
w2_initial = 0.5

# Define a small, identical perturbation for the weights
dw1 = 0.01
dw2 = 0.01
w1_perturbed = w1_initial + dw1
w2_perturbed = w2_initial + dw2

# --- Step 4: Calculate the impact of perturbation for each activation ---
print("This demonstration shows how the activation function's properties affect")
print("the network's sensitivity to weight perturbations.\n")

# Case 1: Standard Tanh (Lipschitz constant approx 1)
L_tanh = np.max(tanh_derivative(np.linspace(-5, 5, 1000)))
output_initial_tanh = simple_nn(x_input, w1_initial, w2_initial, tanh)
output_perturbed_tanh = simple_nn(x_input, w1_perturbed, w2_perturbed, tanh)
output_change_tanh = np.abs(output_perturbed_tanh - output_initial_tanh)

print(f"--- Using standard tanh (Estimated Lipschitz Constant: {L_tanh:.2f}) ---")
print(f"Initial Output Equation: f(x={x_input}; w1={w1_initial}, w2={w2_initial}) = {output_initial_tanh:.5f}")
print(f"Perturbed Output Equation: f(x={x_input}; w1={w1_perturbed:.2f}, w2={w2_perturbed:.2f}) = {output_perturbed_tanh:.5f}")
print(f"Resulting Output Change: |{output_perturbed_tanh:.5f} - {output_initial_tanh:.5f}| = {output_change_tanh:.5f}\n")


# Case 2: Scaled Tanh (Lipschitz constant approx 2)
scaling_factor = 2.0
L_scaled_tanh = np.max(scaled_tanh_derivative(np.linspace(-5, 5, 1000), k=scaling_factor))
# Wrap the scaled_tanh in a lambda to match the expected function signature
activation_scaled_tanh = lambda x: scaled_tanh(x, k=scaling_factor)
output_initial_scaled_tanh = simple_nn(x_input, w1_initial, w2_initial, activation_scaled_tanh)
output_perturbed_scaled_tanh = simple_nn(x_input, w1_perturbed, w2_perturbed, activation_scaled_tanh)
output_change_scaled_tanh = np.abs(output_perturbed_scaled_tanh - output_initial_scaled_tanh)

print(f"--- Using scaled tanh (Estimated Lipschitz Constant: {L_scaled_tanh:.2f}) ---")
print(f"Initial Output Equation: f(x={x_input}; w1={w1_initial}, w2={w2_initial}) = {output_initial_scaled_tanh:.5f}")
print(f"Perturbed Output Equation: f(x={x_input}; w1={w1_perturbed:.2f}, w2={w2_perturbed:.2f}) = {output_perturbed_scaled_tanh:.5f}")
print(f"Resulting Output Change: |{output_perturbed_scaled_tanh:.5f} - {output_initial_scaled_tanh:.5f}| = {output_change_scaled_tanh:.5f}\n")

print("Conclusion: The same weight perturbation causes a larger output change for the")
print("activation with the higher Lipschitz constant, demonstrating its key role in stability.")