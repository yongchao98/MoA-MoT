import sys

def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model for symmetry class D with a given number of replicas.
    """
    # The problem specifies two replicas.
    num_replicas = 2

    # For symmetry class D, the bosonic sector of the sigma-model is
    # parametrized by the symmetric space O(2N, 2N) / (O(2N) x O(2N)),
    # where N is the number of replicas.
    # The dimension of this space is given by the formula p * q,
    # where p = 2N and q = 2N.

    p = 2 * num_replicas
    q = 2 * num_replicas

    # The number of non-Grassman variables is the dimension of this space.
    num_variables = p * q

    print(f"Symmetry Class: D")
    print(f"Number of replicas (N): {num_replicas}")
    print("The manifold for the bosonic variables is O(2N, 2N) / (O(2N) x O(2N)).")
    print("Its dimension is calculated as p * q.")
    print(f"p = 2 * N = 2 * {num_replicas} = {p}")
    print(f"q = 2 * N = 2 * {num_replicas} = {q}")
    print("\nFinal equation for the number of variables:")
    print(f"{p} * {q} = {num_variables}")

if __name__ == '__main__':
    solve_class_d_variables()
