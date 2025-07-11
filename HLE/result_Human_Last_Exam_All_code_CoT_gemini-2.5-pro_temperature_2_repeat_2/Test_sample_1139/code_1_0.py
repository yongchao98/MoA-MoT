import math

def calculate_non_grassman_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model for disordered systems of symmetry class D with n replicas.
    """
    # Number of replicas as specified in the problem
    n = 2

    # In symmetry class D, the target space is the supermanifold
    # OSp(4n|4n) / (OSp(2n|2n) x OSp(2n|2n)).
    # The number of non-Grassman (bosonic) variables is given by the formula: 8 * n^2
    coefficient = 8

    # Calculate the number of variables
    num_variables = coefficient * (n ** 2)

    # Print the explanation and the calculation
    print("The number of non-Grassman variables for the supersymmetric sigma-model")
    print("with n replicas for symmetry class D is given by the formula: 8 * n^2.")
    print(f"\nFor n = {n} replicas, the calculation is:")
    print(f"{coefficient} * {n}**2 = {num_variables}")

if __name__ == '__main__':
    calculate_non_grassman_variables()
