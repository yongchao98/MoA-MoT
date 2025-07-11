import numpy as np

def solve():
    """
    Constructs a set of 9 vectors in C^6 and verifies they meet the specified angle conditions.
    """
    # Set numpy to print complex numbers nicely for readability
    np.set_printoptions(precision=3, suppress=True)

    # 1. Define the 6 standard basis vectors in C^6
    e1 = np.array([1, 0, 0, 0, 0, 0], dtype=complex)
    e2 = np.array([0, 1, 0, 0, 0, 0], dtype=complex)
    e3 = np.array([0, 0, 1, 0, 0, 0], dtype=complex)
    e4 = np.array([0, 0, 0, 1, 0, 0], dtype=complex)
    e5 = np.array([0, 0, 0, 0, 1, 0], dtype=complex)
    e6 = np.array([0, 0, 0, 0, 0, 1], dtype=complex)

    # 2. Define the three additional vectors
    v7 = 0.5 * (e1 + e2 + e3 + e4)
    v8 = 0.5 * (e1 + e2 + e5 + e6)
    v9 = 0.5 * (e3 + e4 + e5 + e6)

    # 3. Create a list of all 9 vectors
    vectors = [e1, e2, e3, e4, e5, e6, v7, v8, v9]
    vector_names = [f"e{i+1}" for i in range(6)] + [f"v{i+7}" for i in range(3)]
    
    print(f"Constructed a set of {len(vectors)} vectors in C^6.")
    print("Verifying conditions...\n")
    
    pi = np.pi
    num_vectors = len(vectors)
    has_orthogonal_pair = False
    
    # 4. Iterate through all unique pairs and check angles
    for i in range(num_vectors):
        for j in range(i + 1, num_vectors):
            v_i = vectors[i]
            v_j = vectors[j]

            # Calculate the Hermitian inner product
            inner_product = np.vdot(v_i, v_j)
            
            # Calculate magnitudes (norms)
            norm_i = np.linalg.norm(v_i)
            norm_j = np.linalg.norm(v_j)
            
            # Calculate the cosine of the angle
            cos_alpha = np.abs(inner_product) / (norm_i * norm_j)
            
            # Convert cosine to angle in radians
            angle = np.arccos(cos_alpha)
            
            # Check if angle is approximately pi/2 or pi/3
            is_pi_over_2 = np.isclose(angle, pi / 2)
            is_pi_over_3 = np.isclose(angle, pi / 3)

            if is_pi_over_2:
                angle_desc = "pi/2"
                has_orthogonal_pair = True
            elif is_pi_over_3:
                angle_desc = "pi/3"
            else:
                # This case should not be reached for our construction
                angle_desc = "INVALID"

            print(f"Angle between {vector_names[i]} and {vector_names[j]}:")
            # The 'final equation' can be interpreted as showing the calculation steps
            print(f"  |(v,w)| / (|v||w|) = {np.abs(inner_product):.3f} / ({norm_i:.3f}*{norm_j:.3f}) = {cos_alpha:.3f}")
            print(f"  cos(alpha) = {cos_alpha:.3f} => alpha is approx. {angle_desc} ({angle*180/pi:.1f} degrees)")
            print("-" * 20)

    if has_orthogonal_pair:
        print("\nVerification successful: All pairs have an angle of pi/2 or pi/3,")
        print("and at least one orthogonal pair exists.")
    else:
        print("\nVerification failed: No orthogonal pair was found.")

    print("\nThe largest possible number of such vectors is 9.")

solve()