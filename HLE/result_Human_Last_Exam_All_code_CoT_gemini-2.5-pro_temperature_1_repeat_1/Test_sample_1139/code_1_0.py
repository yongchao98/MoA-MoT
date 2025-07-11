import math

def main():
    """
    Calculates the number of non-Grassman variables needed to parametrize
    the bosonic sector of the supersymmetric sigma-model with two replicas
    for disordered systems of symmetry class D.
    """
    # Number of replicas as specified in the problem
    N = 2

    # The number of variables corresponds to the dimension of the symmetric space
    # O(2N) / (O(N) x O(N)).
    # We calculate this as dim(O(2N)) - dim(O(N) x O(N)).

    # The dimension of the orthogonal group O(m) is m*(m-1)/2.
    def dim_O(m):
      return m * (m - 1) // 2

    # Parameters for the groups G = O(2N) and K = O(N) x O(N)
    m_G = 2 * N
    m_K = N

    # Calculate the dimension of the total group G = O(4)
    dimension_G = dim_O(m_G)

    # Calculate the dimension of the subgroup K = O(2) x O(2)
    # dim(K) = dim(O(2)) + dim(O(2))
    dimension_K = dim_O(m_K) + dim_O(m_K)

    # The final result is the difference between the dimensions.
    result = dimension_G - dimension_K

    print(f"The calculation is based on the dimension of the symmetric space O({m_G}) / (O({m_K}) x O({m_K})).")
    print(f"The formula for the dimension is: dim(O({m_G})) - (dim(O({m_K})) + dim(O({m_K}))).")
    print("")
    print("Step 1: Calculate the dimension of the total group G = O(4).")
    print(f"dim(O({m_G})) = {m_G} * ({m_G} - 1) / 2 = {dimension_G}")
    print("")
    print("Step 2: Calculate the dimension of the subgroup K = O(2) x O(2).")
    print(f"dim(O({m_K})) = {m_K} * ({m_K} - 1) / 2 = {dim_O(m_K)}")
    print(f"dim(K) = {dim_O(m_K)} + {dim_O(m_K)} = {dimension_K}")
    print("")
    print("Step 3: Calculate the final dimension of the manifold.")
    print(f"Result = dim(O({m_G})) - dim(K)")
    print(f"Result = {dimension_G} - {dimension_K} = {result}")


if __name__ == "__main__":
    main()