def solve_fire_spread_r0f():
    """
    This function defines the variables from the fire spread model
    and prints the derived expression for R0f.
    """
    # Define the symbols for the variables as strings
    b = "b"
    pg = "pg"
    tau = "𝜏"
    c_var = "c"
    pt = "pt"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"

    # Construct the numerator and denominator parts of the expression
    # R0f = (b * pg * 𝜏 * c * pt) / ((𝛾t + 𝜇t) * (𝜏 + 𝜇g) * 𝜇g)
    
    numerator_vars = [b, pg, c_var, pt, tau]
    numerator_str = " * ".join(numerator_vars)
    
    denominator_part1 = f"({gamma_t} + {mu_t})"
    denominator_part2 = f"({tau} + {mu_g})"
    denominator_part3 = f"{mu_g}"
    denominator_str = f"{denominator_part1} * {denominator_part2} * {denominator_part3}"

    # Print the final expression for R0f
    print("The expression for R0f is:")
    print(f"R0f = ( {numerator_str} ) / ( {denominator_str} )")
    
    # As requested, output each variable (symbol) used in the final equation
    print("\nWhere the variables are:")
    print(f"{b}: the burning rate of trees (contacts with grass per day)")
    print(f"{pg}: the probability that a contact between a burning tree and grass ignites the grass")
    print(f"{c_var}: the burning rate of dry grass (contacts with trees per day)")
    print(f"{pt}: the probability that a contact between burning grass and a tree ignites the tree")
    print(f"{tau} (𝜏): the rate at which ignited grass becomes intensely burning")
    print(f"{gamma_t} (𝛾t): the fire extinguishing rate for burning trees")
    print(f"{mu_t} (𝜇t): the rate at which trees naturally die")
    print(f"{mu_g} (𝜇g): the per-day rate at which grass naturally dies")

solve_fire_spread_r0f()