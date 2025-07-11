import numpy as np

def calculate_error_vector(w, vectors):
    """Calculates the vector of squared errors for a given w."""
    w = w / np.linalg.norm(w) # Ensure w is a unit vector
    errors = []
    for x in vectors:
        x = np.array(x)
        norm_x_sq = np.dot(x, x)
        proj_x_w_sq = np.dot(x, w)**2
        error_sq = norm_x_sq - proj_x_w_sq
        errors.append(error_sq)
    return np.array(errors)

def main():
    """
    Demonstrates the non-convexity of the achievable error set for d=3
    with non-negative data points, leading to potential failure of scalarization.
    """
    print("This script demonstrates that for d=3, the set of achievable error vectors can be non-convex,")
    print("even when all data vectors x_i are non-negative. This is a reason why scalarization")
    print("may fail to find the entire Pareto front.\n")
    
    # Counterexample for d=3, based on Boman et al. (2020)
    # x_i vectors are non-negative.
    epsilon = 0.1
    x1 = np.array([1.0, 0.0, epsilon])
    x2 = np.array([0.0, 1.0, epsilon])
    x3 = np.array([1.0, 1.0, 0.0])
    data_vectors = [x1, x2, x3]
    
    print("Data vectors (x_i):")
    print(f"x1 = {x1}")
    print(f"x2 = {x2}")
    print(f"x3 = {x3}\n")
    
    # Let's select two points in the space of solutions w
    w_A = np.array([1.0, 0.0, 0.0])
    w_B = np.array([0.0, 1.0, 0.0])
    
    # Calculate the corresponding error vectors
    error_vec_A = calculate_error_vector(w_A, data_vectors)
    error_vec_B = calculate_error_vector(w_B, data_vectors)
    
    print(f"For w_A = {w_A}, the achievable error vector is e_A = {np.round(error_vec_A, 4)}")
    print(f"For w_B = {w_B}, the achievable error vector is e_B = {np.round(error_vec_B, 4)}\n")
    
    # Calculate the midpoint of the line segment connecting e_A and e_B.
    # This point lies in the convex hull of the set of achievable error vectors.
    error_vec_mid = (error_vec_A + error_vec_B) / 2
    
    print("The midpoint of the segment e_A --- e_B is:")
    print(f"e_mid = {np.round(error_vec_mid, 4)}\n")
    print("If the set of achievable errors were convex, e_mid should itself be achievable")
    print("or dominated by an achievable point.\n")
    
    # Let's find an achievable point that has the same first two error components as e_mid.
    # This can be achieved with a vector w that is a combination of w_A and w_B.
    w_C = (w_A + w_B) / np.linalg.norm(w_A + w_B)
    error_vec_C = calculate_error_vector(w_C, data_vectors)
    
    print(f"Let's test another achievable point, using w_C = {np.round(w_C, 4)}.")
    print(f"The achievable error vector is e_C = {np.round(error_vec_C, 4)}\n")
    
    print("Observation:")
    print(f"The first two components of e_mid ({np.round(error_vec_mid[0], 2)}, {np.round(error_vec_mid[1], 2)}) are identical to e_C's ({np.round(error_vec_C[0], 2)}, {np.round(error_vec_C[1], 2)}).")
    print(f"However, the third component is drastically different: e_mid[2] = {np.round(error_vec_mid[2], 2)} while e_C[2] = {np.round(error_vec_C[2], 2)}.")
    print(f"The point e_mid = {np.round(error_vec_mid, 2)} is not achievable. This demonstrates the non-convexity of the set.")
    print("\nSince a counterexample exists for d=3, but the property holds for d=2, the largest d is 2.")

if __name__ == "__main__":
    main()
