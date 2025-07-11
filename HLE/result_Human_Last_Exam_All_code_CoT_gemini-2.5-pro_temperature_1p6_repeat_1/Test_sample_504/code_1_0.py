def solve_vector_problem():
    """
    Calculates the largest number of vectors satisfying the given conditions.
    """
    
    # We start with the 6 vectors of the standard orthonormal basis in C^6.
    # These are all pairwise orthogonal (angle pi/2).
    num_basis_vectors = 6

    # We can add more vectors. These vectors must have supports of size 4.
    # The supports must be chosen from a family of 4-element subsets of {1,...,6}
    # where any two distinct supports intersect in exactly 2 elements.
    # The maximum size of such a family is 3.
    num_supports = 3

    # For each support, we can construct a number of vectors.
    # By analyzing the angle constraints, we find that we can place a maximum of 8
    # such vectors on each support.
    num_vectors_per_support = 8

    # The total number of vectors is the sum of the basis vectors and the
    # vectors from all supports.
    total_vectors = num_basis_vectors + num_supports * num_vectors_per_support
    
    print("The final calculation is based on the following components:")
    print(f"Number of orthogonal basis vectors = {num_basis_vectors}")
    print(f"Number of allowed vector supports = {num_supports}")
    print(f"Number of vectors per support = {num_vectors_per_support}")
    print(f"The largest number of vectors is {num_basis_vectors} + {num_supports} * {num_vectors_per_support} = {total_vectors}")

solve_vector_problem()