def display_r0f_expression():
    """
    This function prints the derived symbolic expression for R0f.
    It breaks down the formula into its components for clarity.
    """

    # --- Component 1: Number of infectious grass patches from 1 tree ---

    # Average burning period of a tree
    tree_burning_period = "1 / (𝛾t + 𝜇t)"

    # Rate of grass ignition by a tree
    grass_ignition_rate_from_tree = "b * pg"
    
    # Number of grass patches exposed by a single tree
    exposed_grass_from_tree = f"({grass_ignition_rate_from_tree}) * ({tree_burning_period})"
    
    # Probability of a grass patch surviving its latent period to become infectious
    grass_survival_probability = "𝜏 / (𝜏 + 𝜇g)"

    # Total infectious grass patches created by one tree
    infectious_grass_from_tree = f"({exposed_grass_from_tree}) * ({grass_survival_probability})"


    # --- Component 2: Number of trees ignited by 1 infectious grass patch ---
    
    # Average infectious period of a grass patch
    grass_infectious_period = "1 / 𝜇g"
    
    # Rate of tree ignition by a grass patch
    tree_ignition_rate_from_grass = "c * pt"
    
    # Total new trees ignited by one infectious grass patch
    trees_from_grass = f"({tree_ignition_rate_from_grass}) * ({grass_infectious_period})"

    
    # --- Final R0f Expression ---
    
    # R0f is the product of the two main components
    print("The expression for R0f is derived by multiplying the number of infectious grass patches created by one tree")
    print("by the number of new trees ignited by one of those infectious grass patches.")
    print("-" * 50)
    print(f"R0f = (Infectious Grass from 1 Tree) * (Trees from 1 Grass Patch)")
    print(f"R0f = [{infectious_grass_from_tree}] * [{trees_from_grass}]")
    print("-" * 50)

    # Simplified final formula
    # Numerator contains all the individual rate and probability terms
    numerator = "b * pg * c * pt * 𝜏"
    # Denominator contains all the rate sums and lifetime terms
    denominator = "(𝛾t + 𝜇t) * (𝜏 + 𝜇g) * 𝜇g"
    
    print("The final simplified expression for R0f is:")
    print(f"      {numerator}")
    print(f"R0f = {'-' * (len(numerator) + 6)}")
    print(f"      {denominator}")

if __name__ == '__main__':
    display_r0f_expression()
