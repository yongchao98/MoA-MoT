def show_alpha_scaling_relation():
    """
    This function prints the quantitative scaling relationship between the specific
    heat critical exponent (α) and the spatial dimensionality (d) for d < 4.

    The relationship is derived from the epsilon expansion in Renormalization
    Group theory for an O(n)-symmetric scalar field theory.
    """

    # Symbolic representations
    alpha_sym = "α"
    d_sym = "d"
    n_sym = "n"  # Number of components of the scalar field

    # Constants from the first-order epsilon expansion
    const_four = 4
    const_two = 2
    const_eight = 8

    # Introduction explaining the context
    print("In scalar field theory, the scaling of the specific heat exponent (α)")
    print("with spatial dimensionality (d) below the upper critical dimension (d=4)")
    print("is determined using the epsilon expansion (ε = 4 - d).")
    print("\nFor a general O(n)-symmetric model, the first-order relationship is:")

    # Print the final equation with all numbers shown explicitly
    equation = (
        f"{alpha_sym} ≈ ( ({const_four} - {n_sym}) / ({const_two} * ({n_sym} + {const_eight})) ) * ({const_four} - {d_sym})"
    )
    print("\n" + "="*len(equation))
    print(equation)
    print("="*len(equation) + "\n")

    print(f"Where:")
    print(f"  {alpha_sym} is the critical exponent for specific heat.")
    print(f"  {d_sym} is the spatial dimensionality.")
    print(f"  {n_sym} is the number of components of the order parameter (e.g., n=1 for the Ising model).")

# Execute the function to display the result
show_alpha_scaling_relation()
