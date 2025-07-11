import torch
import torch.nn as nn
import math

def demonstrate_initialization(fan_in, fan_out, magnitude_scale_factor):
    """
    Demonstrates initializing a weight tensor with a specific magnitude.
    
    In theoretical analyses, the variance of weights is often scaled.
    - Standard He/Kaiming init: variance = 2 / fan_in
    - Standard Glorot/Xavier init: variance = 2 / (fan_in + fan_out)
    
    The 'magnitude_scale_factor' allows us to deviate from these standard
    initializations, which is the key parameter for switching between
    theoretical regimes (e.g., lazy vs. rich training).
    """
    
    # Using Kaiming Uniform as a base
    std = math.sqrt(2.0 / fan_in)
    
    # The scale factor directly controls the magnitude of the initial weights
    # A small factor might lead to a "lazy" regime, a large one to a "rich" regime.
    scaled_std = std * magnitude_scale_factor
    
    # Calculate the bounds for a uniform distribution based on the standard deviation
    # For a uniform distribution U(-a, a), the variance is a^2 / 3. So, a = sqrt(3 * variance).
    # a = sqrt(3) * std
    bound = math.sqrt(3.0) * scaled_std
    
    # Create a weight tensor
    weights = torch.zeros(fan_out, fan_in)
    # Initialize from a uniform distribution with the calculated bound
    nn.init.uniform_(weights, -bound, bound)
    
    actual_std = torch.std(weights).item()
    
    print(f"--- Initialization with scale factor: {magnitude_scale_factor} ---")
    print(f"Network property being set: The magnitude of weight initialization.")
    print(f"This is controlled by a scale factor applied to a standard initialization scheme.")
    print(f"The equation for the bound of the uniform distribution is: bound = sqrt(3 * 2 / fan_in) * scale_factor")
    print(f"For fan_in={fan_in} and scale_factor={magnitude_scale_factor}, the calculated bound is: {bound:.4f}")
    print(f"The first few initialized weights are:\n{weights[:2, :4]}")
    print(f"The standard deviation of the initialized weights is: {actual_std:.4f}")
    print("-" * 20)

# Example parameters for a layer in a neural network
input_features = 512
output_features = 512

# Small magnitude (closer to NTK/lazy regime)
demonstrate_initialization(input_features, output_features, magnitude_scale_factor=1.0)

# Large magnitude (closer to mean-field/rich regime)
# This requires a larger scale factor. For instance, sqrt(fan_in) would be a large factor.
large_scale_factor = math.sqrt(input_features / 2) # Factor to make variance O(1) instead of O(1/fan_in)
demonstrate_initialization(input_features, output_features, magnitude_scale_factor=large_scale_factor)
