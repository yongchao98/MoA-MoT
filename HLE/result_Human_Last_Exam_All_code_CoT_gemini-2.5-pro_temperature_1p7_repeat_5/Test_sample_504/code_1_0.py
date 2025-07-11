import numpy as np
import math

# Plan:
# 1. Construct the 21 vectors as normalized differences of standard basis vectors in R^7.
#    These vectors lie in a 6-dimensional subspace where the sum of coordinates is zero.
# 2. Compute the Gram matrix (matrix of all pairwise inner products).
# 3. Check if the absolute values of off-diagonal elements are either 0 or 0.5.
# 4. Verify that at least one off-diagonal element is 0.
# 5. Print the final result based on the verification and the theoretical bound.

def main():
    """
    Finds the largest number of vectors in C^6 under the given angle conditions.
    """
    # Parameters from the theoretical bound
    d = 6  # dimension of the space
    # The set of squared inner products between distinct vectors is {0, (1/2)^2} = {0, 0.25}
    m = 2  # number of allowed squared inner products
    
    # We construct vectors in a 7D ambient space that lie in a 6D subspace.
    ambient_dim = 7
    num_vectors = math.comb(ambient_dim, 2)

    # Construct the vectors: (e_i - e_j)
    vectors = []
    for i in range(ambient_dim):
        for j in range(i + 1, ambient_dim):
            v = np.zeros(ambient_dim)
            v[i] = 1
            v[j] = -1
            vectors.append(v)
    vectors = np.array(vectors)

    # Normalize the vectors
    unit_vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]

    # Calculate the Gram matrix of inner products
    gram_matrix = np.dot(unit_vectors, unit_vectors.T)
    # Round to handle floating point inaccuracies
    gram_matrix = np.round(gram_matrix, decimals=10)

    # Check the conditions
    angle_conditions_met = True
    orthogonal_pair_found = False
    
    for i in range(num_vectors):
        for j in range(i + 1, num_vectors):
            abs_inner_product = abs(gram_matrix[i, j])
            
            # Check angle condition: |(v,w)| should be 0 or 0.5
            if not math.isclose(abs_inner_product, 0.0) and not math.isclose(abs_inner_product, 0.5):
                angle_conditions_met = False
                break
            
            # Check for an orthogonal pair
            if math.isclose(abs_inner_product, 0.0):
                orthogonal_pair_found = True
        if not angle_conditions_met:
            break

    print(f"Number of vectors constructed: {num_vectors}")
    
    if angle_conditions_met:
        print("The angle conditions (angles are pi/2 or pi/3) are satisfied for all pairs.")
    else:
        print("The angle conditions failed.")

    if orthogonal_pair_found:
        print("At least one pair of orthogonal vectors was found.")
    else:
        print("No orthogonal pair was found.")

    if angle_conditions_met and orthogonal_pair_found:
        print("\nThe construction of 21 vectors is valid.")
        print("This matches the theoretical upper bound k <= C(d+m-1, m).")
        upper_bound = math.comb(d + m - 1, m)
        print(f"k <= C({d}+{m}-1, {m}) = C({d+m-1}, {m}) = (({d+m-1} * {d+m-2}) / 2) = ({7} * {6}) / 2 = {upper_bound}")
        print("\nThus, the largest number of such vectors is 21.")

if __name__ == '__main__':
    main()
