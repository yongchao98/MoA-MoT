import numpy as np
from numpy.linalg import solve, lstsq

def main():
    """
    Analyzes LIME feature importances for a given model and two different scenarios.
    """
    # Define the model function
    lookup_table = {1.0: 1.0, 0.0: 0.0}
    def f(input1, input2):
        """The black-box model to explain."""
        return lookup_table.get(input1, input1 * 0 + input2 * 0.5 + 0.5)

    print("Analyzing LIME feature importance for the given model.")
    print("=" * 60)

    # --- Case i) E belongs to the baseline dataset, e.g., E = (0.0, 0.0) ---
    print("Case i): Explaining instance E = (0.0, 0.0)")
    print("LIME creates binary features (z') for whether a feature was perturbed.")
    print("z'_1=1 if input1 not perturbed, 0 otherwise.")
    print("z'_2=1 if input2 not perturbed, 0 otherwise.")
    print("Fitting model: f(z) ~ bias + w1*z'_1 + w2*z'_2\n")

    # We create a small dataset representing the 4 perturbation scenarios from E=(0,0).
    # The baseline value used for replacement is 1.0.
    # z: perturbed instance, z_prime: LIME's binary representation
    # 1. No perturbation:   z = (0,0), z_prime = (1,1). f(z) = 0.0
    # 2. Perturb input1:     z = (1,0), z_prime = (0,1). f(z) = 1.0
    # 3. Perturb input2:     z = (0,1), z_prime = (1,0). f(z) = 0.0
    # 4. Perturb both:       z = (1,1), z_prime = (0,0). f(z) = 1.0
    
    # Setup the linear system A * x = y, for x = [bias, w1, w2]^T
    A_i = np.array([
        # bias, z'_1, z'_2
        [1,    1,      1],     # Corresponds to z_prime=(1,1)
        [1,    0,      1],     # Corresponds to z_prime=(0,1)
        [1,    1,      0],     # Corresponds to z_prime=(1,0)
        [1,    0,      0]      # Corresponds to z_prime=(0,0)
    ])
    y_i = np.array(
        [0.0, 1.0, 0.0, 1.0]
    )

    coeffs_i = lstsq(A_i, y_i, rcond=None)[0]
    bias_i, w1_i, w2_i = coeffs_i

    print("The final equation for the local linear model is:")
    print(f"f(z) = {bias_i:.2f} + ({w1_i:.2f} * z'_1) + ({w2_i:.2f} * z'_2)")
    
    print(f"\nFeature importances (absolute weights): |w1| = {abs(w1_i):.2f}, |w2| = {abs(w2_i):.2f}")
    if abs(w1_i) > abs(w2_i):
        print("Conclusion for i): input1 is more important.")
    else:
        print("Conclusion for i): input2 is more important.")

    print("\n" + "=" * 60 + "\n")

    # --- Case ii) E does not belong to the baseline, e.g., E = (-1.0, -1.0) ---
    print("Case ii): Explaining instance E = (-1.0, -1.0)")
    print("LIME fits a weighted linear model: f(z) ~ bias + w1*z_1 + w2*z_2")
    print("We model this by solving for the closest perturbations.\n")

    # The prediction for E=(-1,-1) is f(-1,-1) = 0.5*(-1) + 0.5 = 0.0.
    # We simulate LIME's weighted regression by using E and its two closest
    # single-feature perturbations (using baseline value 0.0 for replacement).
    # 1. The instance itself: E = (-1,-1). f(-1,-1) = 0.0
    # 2. Perturb input1 only: z1 = (0, -1). f(0,-1) = 0.0
    # 3. Perturb input2 only: z2 = (-1, 0). f(-1,0) = 0.5
    
    # Setup the linear system A * x = y for x = [bias, w1, w2]^T
    # bias + w1*(-1) + w2*(-1) = 0.0
    # bias + w1*(0)  + w2*(-1) = 0.0
    # bias + w1*(-1) + w2*(0)  = 0.5
    A_ii = np.array([
        # bias, z_1,  z_2
        [1,    -1,   -1],
        [1,     0,   -1],
        [1,    -1,    0]
    ])
    y_ii = np.array([0.0, 0.0, 0.5])

    coeffs_ii = solve(A_ii, y_ii)
    bias_ii, w1_ii, w2_ii = coeffs_ii

    print("The final equation for the local linear model is:")
    print(f"f(z) = {bias_ii:.2f} + ({w1_ii:.2f} * z_1) + ({w2_ii:.2f} * z_2)")

    print(f"\nFeature importances (absolute weights): |w1| = {abs(w1_ii):.2f}, |w2| = {abs(w2_ii):.2f}")
    if abs(w1_ii) > abs(w2_ii):
        print("Conclusion for ii): input1 is more important.")
    else:
        print("Conclusion for ii): input2 is more important.")

if __name__ == '__main__':
    main()
<<<A>>>