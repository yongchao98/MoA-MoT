import math

def main():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model in class D with two replicas.
    """

    # 1. Define constants and the helper function for group dimension
    n_replicas = 2

    def dim_orthogonal_group(m):
      """Calculates the dimension of the orthogonal group O(m)."""
      if m <= 0:
        return 0
      # The dimension of O(m) is m*(m-1)/2
      return m * (m - 1) // 2

    # 2. Main calculation
    # Parameters for the G/K symmetric spaces based on n_replicas
    N_large = 2 * n_replicas
    N_small = n_replicas

    # Calculate dimensions of the relevant Lie groups
    # Numerator group G = O(2n)
    dim_O_2n = dim_orthogonal_group(N_large)
    # Factor of the denominator group K = O(n)
    dim_O_n = dim_orthogonal_group(N_small)
    # Full denominator group K = O(n) x O(n)
    dim_K_factor_group = 2 * dim_O_n

    # The space of bosonic variables decomposes into two identical sectors
    # (M_B and M_F), each of which is a Grassmannian O(2n)/(O(n)xO(n)).
    # We first calculate the dimension of one such sector.
    dim_sector = dim_O_2n - dim_K_factor_group

    # The total number of bosonic variables is the sum of the dimensions
    # of the two sectors.
    total_dim = 2 * dim_sector

    # 3. Print the explanation and step-by-step result
    print("The number of bosonic variables is the dimension of the target manifold's bosonic part.")
    print("For class D, this space decomposes into two identical sectors M_B and M_F.")
    print(f"Each sector is the Grassmannian manifold O(2n)/(O(n)xO(n)), with n={n_replicas} replicas.")
    print("-" * 30)
    print(f"Step 1: Calculate the dimension of the numerator group O(2n) = O({N_large}).")
    print(f"   dim(O({N_large})) = {N_large}*({N_large}-1)/2 = {dim_O_2n}")

    print(f"\nStep 2: Calculate the dimension of the denominator group O(n)xO(n) = O({N_small})xO({N_small}).")
    print(f"   dim(O({N_small})) = {N_small}*({N_small}-1)/2 = {dim_O_n}")
    print(f"   dim(O({N_small})xO({N_small})) = dim(O({N_small})) + dim(O({N_small})) = {dim_O_n} + {dim_O_n} = {dim_K_factor_group}")

    print(f"\nStep 3: Calculate the dimension of a single sector.")
    print(f"   dim(Sector) = dim(O({N_large})) - dim(O({N_small})xO({N_small})) = {dim_O_2n} - {dim_K_factor_group} = {dim_sector}")

    print("\nStep 4: The total number of variables is dim(M_B) + dim(M_F).")
    print("The final result is given by the equation:")
    print(f"{dim_sector} + {dim_sector} = {total_dim}")


if __name__ == "__main__":
    main()