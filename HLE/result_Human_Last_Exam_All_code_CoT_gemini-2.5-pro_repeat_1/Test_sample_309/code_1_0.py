def calculate_max_marginalized_landmarks():
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem while ensuring solvability.
    """
    # For a bundle adjustment problem to be well-posed, we generally need
    # at least 3 landmarks (N>=3) and 2 cameras (M>=2).
    # Let's use some example values that satisfy these conditions.
    N = 20  # Total number of landmarks
    M = 10  # Total number of cameras

    print(f"Given a bundle adjustment problem with N = {N} landmarks and M = {M} cameras.")
    print("We want to find the maximum number of landmarks that can be marginalized using the Schur complement.")
    print("\nReasoning:")
    print("The solvability of the system depends on having a stable geometric reference.")
    print("A minimum of 3 non-collinear landmarks are required to form a stable 3D frame and uniquely determine camera poses.")
    print("When we marginalize 'k' landmarks, we are left solving a reduced system with 'N - k' landmarks.")
    print("To ensure this reduced system is solvable, the number of landmarks kept must be at least 3.")
    print("Therefore, the condition is: N - k >= 3")
    print("Solving for the maximum k gives: k <= N - 3")

    # The maximum number of landmarks that can be marginalized is N - 3.
    max_k = N - 3

    print("\nResulting Equation:")
    # The final output prints each number in the equation as requested.
    print(f"Max marginalized landmarks = {N} - 3 = {max_k}")


# Execute the function
calculate_max_marginalized_landmarks()