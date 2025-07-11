import math

def main():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model in symmetry class D with a specified number of replicas.
    """
    # The problem specifies two replicas for the disordered system.
    N = 2

    # The formula for the dimension of the orthogonal group O(n) is n*(n-1)/2.
    def dim_O(n):
        return n * (n - 1) // 2

    # The bosonic manifold for symmetry class D is O(2N) / (O(N) x O(N)).
    # We calculate the dimension of the numerator group, O(2N).
    numerator_n = 2 * N
    dim_numerator = dim_O(numerator_n)

    # We calculate the dimension of the denominator group, O(N) x O(N).
    # Its dimension is dim(O(N)) + dim(O(N)).
    denominator_n = N
    dim_denominator_single = dim_O(denominator_n)
    dim_denominator_total = 2 * dim_denominator_single

    # The number of variables is the dimension of the manifold, which is
    # the difference between the dimensions of the numerator and denominator groups.
    num_variables = dim_numerator - dim_denominator_total

    print("Calculation for the number of bosonic variables in the sigma-model for class D with N=2 replicas:")
    print("-" * 80)
    print(f"The target manifold is the symmetric space O(2N) / (O(N) x O(N)), with N = {N}.")
    print(f"The dimension is calculated as: dim(O(2N)) - (dim(O(N)) + dim(O(N))).")
    print("-" * 80)

    print(f"1. Dimension of the numerator group O({numerator_n}):")
    print(f"   dim(O({numerator_n})) = {numerator_n} * ({numerator_n} - 1) / 2 = {dim_numerator}")

    print(f"\n2. Dimension of the denominator group O({denominator_n}) x O({denominator_n}):")
    print(f"   dim(O({denominator_n})) = {denominator_n} * ({denominator_n} - 1) / 2 = {dim_denominator_single}")
    print(f"   Total dimension = {dim_denominator_single} + {dim_denominator_single} = {dim_denominator_total}")

    print("-" * 80)
    print("3. Final Result:")
    print("The number of non-Grassman variables is the difference in dimensions.")
    # The final equation with each number printed explicitly.
    print(f"Number of variables = {dim_numerator} - {dim_denominator_total} = {num_variables}")

if __name__ == "__main__":
    main()
