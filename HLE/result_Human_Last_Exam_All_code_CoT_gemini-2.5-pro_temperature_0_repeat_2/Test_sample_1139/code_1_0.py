import math

def calculate_variables_class_d(n_replicas):
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model of symmetry class D.

    The bosonic sector's target manifold for class D with n replicas is O(2n)/U(n).
    The number of variables is the dimension of this manifold.
    """
    print(f"For symmetry class D with n={n_replicas} replicas, the target manifold for the bosonic sector is O(2*n) / U(n).")
    
    # Parameters for the groups
    k_O = 2 * n_replicas
    k_U = n_replicas
    
    print(f"For n={n_replicas}, this becomes the manifold O({k_O}) / U({k_U}).")
    print("\nThe number of variables is the dimension of this manifold, calculated as dim(O(k)) - dim(U(k)).")

    # Dimension of O(k) = k * (k - 1) / 2
    dim_O = (k_O * (k_O - 1)) // 2
    print(f"\nFirst, we calculate the dimension of O({k_O}):")
    print(f"dim(O({k_O})) = ({k_O} * ({k_O} - 1)) / 2 = {dim_O}")

    # Dimension of U(k) = k^2
    dim_U = k_U ** 2
    print(f"\nNext, we calculate the dimension of U({k_U}):")
    print(f"dim(U({k_U})) = {k_U}^2 = {dim_U}")

    # Dimension of the quotient space
    num_variables = dim_O - dim_U
    
    print("\nFinally, we find the difference to get the number of non-Grassman variables:")
    print(f"Number of variables = dim(O({k_O})) - dim(U({k_U}))")
    print(f"The final equation is: {dim_O} - {dim_U} = {num_variables}")

# Number of replicas given in the problem
n = 2
calculate_variables_class_d(n)