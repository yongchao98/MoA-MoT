import numpy as np

def print_equation(S_new, E_out2, LHS, E_out1):
    """Helper function to print the final equation."""
    print("--- Testing if the theory S_new * E_out2 = E_out1 holds ---\n")
    
    print("Matrix S_new (B @ T):")
    print(np.round(S_new, 2))
    print("\n* (multiplied by)\n")
    
    print("Vector E_out2 (calculated using old inverse):")
    print(np.round(E_out2, 2))
    print("\n= (gives result, LHS)\n")
    
    print("Result Vector (LHS):")
    print(np.round(LHS, 2))
    print("\n----------------------------------------------------")
    
    print("Original Output Vector (RHS, E_out1):")
    print(np.round(E_out1, 2))
    print("\n")

    # Final check and conclusion
    if np.allclose(LHS, E_out1):
        print("Conclusion: The equation holds. The theory is still valid.")
    else:
        print("Conclusion: The equation does not hold. The theory fails.")
        print("As you can see, the Result Vector (LHS) is not equal to the Original Output Vector (RHS).")


# Set up a random number generator for reproducibility
rng = np.random.default_rng(0)

# 1. Define the optical components as 2x2 Jones Matrices.
# We are ignoring spatial effects to focus on polarization, which is where the birefringent medium acts.

# A random medium `T` that couples polarizations.
# We create a random complex matrix to simulate this.
T = rng.random((2, 2)) + 1j * rng.random((2, 2))
T_inv = np.linalg.inv(T)

# A birefringent medium `B`, simulated as a quarter-wave plate with its fast axis at 0 degrees.
# This introduces a phase shift between horizontal and vertical components.
phi = np.pi / 2  # 90-degree phase shift
B = np.array([[1, 0], [0, np.exp(1j * phi)]])

# 2. Define the input beam's polarization (e.g., horizontal).
E_in = np.array([1, 0])

# 3. The new system includes the birefringent plate. Let's assume light passes through T then B.
S_new = B @ T

# 4. Propagate the input beam through the NEW system to get Output 1.
# This is the reference output we are trying to achieve.
E_out1 = S_new @ E_in

# 5. Now, apply the procedure from the original theory to find the "necessary input".
# The procedure uses the inverse of the ORIGINAL random medium (T_inv),
# not the inverse of the new, complete system (S_new_inv).
# This is the crucial step where the error is introduced.
E_out2 = T_inv @ E_out1

# 6. According to the theory, propagating E_out2 through the new system S_new
# should result in E_out1. Let's test this.
# LHS (Left-Hand Side of equation) = S_new * E_out2
LHS = S_new @ E_out2
# RHS (Right-Hand Side of equation) = E_out1

# 7. Print the components of the final equation to show the result.
print_equation(S_new, E_out2, LHS, E_out1)
