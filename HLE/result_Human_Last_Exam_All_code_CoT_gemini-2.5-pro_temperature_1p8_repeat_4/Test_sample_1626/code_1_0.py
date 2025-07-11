import numpy as np
import matplotlib.pyplot as plt

def demonstrate_tanh_stretching(N=20, k=2.0):
    """
    Demonstrates how the tanh function is used for grid stretching.

    Args:
        N (int): Number of grid points.
        k (float): Stretching factor. Higher k means more points clustered at the ends.
    """
    # 1. Create a uniform grid in "computational space" from -1 to 1
    xi = np.linspace(-1, 1, N)

    # 2. Apply the tanh mapping to get the non-uniform grid in "physical space"
    # The division by tanh(k) ensures the domain remains from -1 to 1
    x = np.tanh(k * xi) / np.tanh(k)

    # 3. Print the resulting physical grid points
    print(f"Demonstration of a grid with {N} points and a stretching factor k={k}:")
    print("Uniform computational grid points (xi):")
    print(np.round(xi, 4))
    print("\nStretched physical grid points (x) using tanh:")
    print(np.round(x, 4))

    # 4. Plot for visualization
    plt.figure(figsize=(10, 2))
    plt.plot(x, np.zeros_like(x), 'o-', label='Stretched Grid (tanh)')
    plt.plot(xi, np.zeros_like(xi) - 0.1, '|', markersize=20, label='Uniform Grid')
    plt.yticks([])
    plt.title(f'Grid Stretching using tanh (k={k})')
    plt.xlabel('Physical Coordinate (x)')
    plt.legend(loc='lower right')
    plt.ylim(-0.5, 0.5)
    plt.grid(True, axis='x', linestyle=':')
    # plt.show() # Disabled to prevent popping up a window in a server environment
    print("\nA plot would visually show points clustered near -1 and 1.")


# The core question is about the *name* of the function.
# The code above demonstrates *why* tanh is the answer.
# Let's provide the answer as requested.

answer_function_name = "tanh"

print("\n---")
print("The mathematical function that governs grid spacing for resolution consistency from the given choices is:")
print(answer_function_name)
