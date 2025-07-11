def define_set_M():
    """
    This function prints the definition of the set M required to prove the
    existence and uniqueness of the solution to the given boundary value problem
    using the Banach fixed-point theorem.
    """

    # Define the components of the set M's definition
    space = "C([0, 1])"
    description = "the space of continuous functions on the interval [0, 1]"
    
    # The numbers used in the definition
    num_0 = 0
    num_1 = 1

    # Define the conditions for a function u to be in M
    boundary_cond_1 = f"u({num_0}) = {num_0}"
    boundary_cond_2 = f"u({num_1}) = {num_0}"
    inequality_cond = f"u(x) <= {num_0} for all x in [{num_0}, {num_1}]"

    # Construct the final definition string
    set_definition = (
        f"M = {{u ∈ {space} | {boundary_cond_1}, {boundary_cond_2}, and {inequality_cond}}}"
    )

    print("To prove the existence and uniqueness of a global solution using the Banach fixed-point theorem,")
    print("the appropriate set M to define is:")
    print("\n" + "="*80)
    print(set_definition)
    print("="*80 + "\n")
    print(f"where {space} is {description}.")
    print(f"The numbers used in this definition are {num_0} and {num_1}.")

# Execute the function to print the answer
define_set_M()