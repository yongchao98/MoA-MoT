def calculate_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model of symmetry class D with a given number of replicas.
    """
    # Number of replicas as specified in the problem
    N = 2

    # In the supersymmetric sigma-model formalism, the bosonic sector for
    # symmetry class D with N replicas is described by the orthogonal group O(2N).
    # The number of non-Grassman variables is the dimension of this manifold.
    print(f"The problem is for symmetry class D with N = {N} replicas.")
    print("The manifold for the bosonic sector is the orthogonal group O(2N).")

    # The order of the orthogonal group O(n) is n.
    # For our case, n = 2 * N.
    n = 2 * N
    print(f"For N = {N}, the specific group is O({n}).")
    print("")

    # The dimension of the orthogonal group O(n) is given by the formula:
    # dim(O(n)) = n * (n - 1) / 2
    # This formula counts the number of independent generators of the group.
    print("The number of variables is the dimension of this group.")
    print(f"The dimension of O(n) is calculated using the formula: n * (n - 1) / 2.")
    print(f"Substituting n = {n}, the calculation is:")

    # Perform the calculation step-by-step
    numerator = n * (n - 1)
    dimension = numerator // 2

    print(f"  Dimension = {n} * ({n} - 1) / 2")
    print(f"  Dimension = {n} * {n - 1} / 2")
    print(f"  Dimension = {numerator} / 2")
    print(f"  Dimension = {dimension}")
    print("")

    print(f"The final answer is {dimension}.")

if __name__ == "__main__":
    calculate_variables()