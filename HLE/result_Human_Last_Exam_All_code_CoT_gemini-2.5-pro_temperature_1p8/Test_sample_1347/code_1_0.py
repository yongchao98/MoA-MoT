def solve_R0f_expression():
    """
    This function programmatically derives and prints the expression for R0f.
    """

    # Define symbolic representations of the variables
    b = "b"
    c = "c"
    pg = "pg"
    pt = "pt"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"

    # Step 1: Rate of removal for a burning tree (extinguishing + natural death)
    # The average burning period for a single tree is the inverse of this rate.
    tree_burning_period = f"1 / ({gamma_t} + {mu_t})"

    # Step 2: Total number of grass patches ignited by one tree over its lifetime.
    # This is (rate of grass ignition per tree) * (average burning period of a tree)
    grass_ignited_by_tree = f"({b} * {pg}) / ({gamma_t} + {mu_t})"

    # Step 3: Rate of removal for a burning grass patch (natural death only)
    # The average burning period for a grass patch is the inverse of this rate.
    grass_burning_period = f"1 / {mu_g}"
    
    # Step 4: Total number of trees ignited by a single grass patch over its lifetime.
    # This is (rate of tree ignition per grass patch) * (average burning period of a grass patch)
    trees_ignited_by_grass = f"({c} * {pt}) / {mu_g}"
    
    # Step 5: R0f is the product of the two intermediate steps.
    # (Number of grass ignited per tree) * (Number of trees ignited per grass)
    numerator = f"{b} * {c} * {pg} * {pt}"
    denominator = f"({gamma_t} + {mu_t}) * {mu_g}"

    # --- Output Section ---
    print("Derivation of R0f:")
    print("1. A single tree burns for an average period of: " + tree_burning_period)
    print("2. In this time, it ignites this many grass patches: " + grass_ignited_by_tree)
    print("3. Each grass patch then burns for an average period of: " + grass_burning_period)
    print("4. Each burning grass patch ignites this many trees: " + trees_ignited_by_grass)
    
    print("\nMultiplying (2) and (4) gives the final expression for R0f.")
    print("\nFinal Equation:")
    # Printing each variable in the final equation
    print(f"R0f = ({b} * {c} * {pg} * {pt}) / (({gamma_t} + {mu_t}) * {mu_g})")

solve_R0f_expression()