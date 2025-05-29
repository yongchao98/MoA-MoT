import numpy as np

# Define the DenseLayer and LeakyReLU classes
class DenseLayer:
    def __init__(self, num_units, input_shape=None):
        self.num_units = num_units
        self.in_shape = (input_shape, 1)
        self.W = np.random.random((self.in_shape[0], self.num_units))
        self.bias = np.zeros((1, self.num_units))
    
    def forward_pass(self, x):
        return np.dot(x, self.W) + self.bias

class LeakyReLU:
    def __init__(self, alpha=0.01):
        self.alpha = alpha
    
    def forward_pass(self, x):
        return np.where(x < 0, self.alpha * x, x)

# Main function to simulate the process
def main_solution(num_units, input_shape, alpha):
    # Initialize the layers
    dense_layer = DenseLayer(num_units, input_shape)
    leaky_relu_layer = LeakyReLU(alpha)
    
    # Generate a random input
    x = np.random.random((1, input_shape))
    
    # Perform forward pass through the layers
    dense_output = dense_layer.forward_pass(x)
    final_output = leaky_relu_layer.forward_pass(dense_output)
    
    return final_output.tolist()

# Execute the function with assumed parameters
output = main_solution(25, 25, 0.01)
print(output)