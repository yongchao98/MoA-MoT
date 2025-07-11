import numpy as np

def solve_vector_beam_problem():
    """
    Demonstrates why an arbitrary vector beam cannot be generated from a
    scalar input field passing through a fixed linear system.
    """
    # 1. Define the parameters for our discretized system model.
    # We model the beam's cross-section with N pixels.
    N = 4
    # Use a fixed seed for reproducibility of the "random" system.
    np.random.seed(0)

    # 2. Create a random linear system L.
    # The system transforms a single scalar input field 'f' (N-dimensional vector)
    # into a vector output field with two components, g_x and g_y.
    # g_x = L_xx * f
    # g_y = L_yx * f
    # L_xx and L_yx are NxN matrices representing the system's properties.
    L_xx = np.random.rand(N, N) + 1j * np.random.rand(N, N)
    L_yx = np.random.rand(N, N) + 1j * np.random.rand(N, N)

    # The system must be able to resolve the input, so L_xx must be invertible.
    # A random matrix is almost surely invertible.
    if np.linalg.matrix_rank(L_xx) < N:
        print("Error: The system sub-matrix L_xx is singular. Cannot solve.")
        return

    # 3. Define an ARBITRARY target vector beam we want to generate.
    # Both components are chosen randomly to represent an arbitrary target.
    g_x_target = np.random.rand(N) + 1j * np.random.rand(N)
    g_y_target = np.random.rand(N) + 1j * np.random.rand(N)

    print("--- STEP 1: Define an Arbitrary Target Vector Beam ---")
    print(f"Target x-component (g_x_target):\n{np.round(g_x_target, 2)}\n")
    print(f"Target y-component (g_y_target):\n{np.round(g_y_target, 2)}\n")
    print("-" * 50)

    # 4. Find the required input `f` to produce the TARGET x-component.
    # We solve the equation: L_xx * f = g_x_target
    # The solution is: f = inv(L_xx) * g_x_target
    L_xx_inv = np.linalg.inv(L_xx)
    f_required = np.dot(L_xx_inv, g_x_target)

    print("--- STEP 2: Calculate Input 'f' to Achieve Target X-Component ---")
    print("From g_x_target = L_xx * f, we find the required input f.")
    print(f"Required scalar input field (f):\n{np.round(f_required, 2)}\n")
    print("-" * 50)

    # 5. Calculate the y-component that this specific `f` ACTUALLY produces.
    # We use the second system equation: g_y_result = L_yx * f
    g_y_result = np.dot(L_yx, f_required)
    
    # This is the "final equation" we need to evaluate.
    # g_y_result = L_yx * inv(L_xx) * g_x_target

    print("--- STEP 3: Calculate the Resulting Y-Component ---")
    print("This required input 'f' will simultaneously generate a y-component given by g_y = L_yx * f.")
    print("Resulting y-component (g_y_result):")
    # Output each number in the final equation's result vector
    print(np.round(g_y_result, 2))
    print("\n")
    print("-" * 50)

    # 6. Compare the resulting y-component with the target y-component.
    error = np.linalg.norm(g_y_target - g_y_result)
    
    print("--- STEP 4: Compare Result with Target ---")
    print("Is the resulting y-component equal to our arbitrary target y-component?\n")
    print(f"Target y-component:   {np.round(g_y_target, 2)}")
    print(f"Resulting y-component:  {np.round(g_y_result, 2)}")
    print(f"\nDifference (norm) between target and result: {error:.4f}")
    
    print("\nCONCLUSION:")
    print("As the norm of the difference is not zero, the resulting y-component does not match the target.")
    print("This demonstrates that while we can achieve any desired x-component by tailoring the input,")
    print("we are not free to choose the corresponding y-component. Therefore, an arbitrary vector beam")
    print("cannot be generated.")

solve_vector_beam_problem()