def solve_fire_spread_expression():
    """
    This function derives and prints the symbolic expression for R0f,
    the basic reproduction number for fire spreading from trees.
    """
    
    # Symbolic representations of the variables defined in the problem.
    b = "b"
    pg = "pg"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    tau = "𝜏"
    mu_g = "𝜇g"
    c = "c"
    pt = "pt"

    # The expression for R0f is the product of three terms:
    # 1. The number of grass patches ignited by one tree over its lifetime.
    # 2. The probability that an ignited grass patch becomes a spreading fire.
    # 3. The number of trees ignited by one spreading grass fire over its lifetime.
    
    # Term 1: Number of grass fires per tree fire
    # Rate of grass ignition * Lifetime of burning tree
    # (b * pg) * (1 / (𝛾t + 𝜇t))
    term1_numerator = f"{b} * {pg}"
    term1_denominator = f"({gamma_t} + {mu_t})"
    
    # Term 2: Probability of an ignited grass patch becoming a spreading fire
    # This is a competing risks model: becoming spreading (rate 𝜏) vs dying (rate 𝜇g)
    # 𝜏 / (𝜏 + 𝜇g)
    term2_numerator = f"{tau}"
    term2_denominator = f"({tau} + {mu_g})"
    
    # Term 3: Number of tree fires per spreading grass fire
    # Rate of tree ignition * Lifetime of spreading grass
    # (c * pt) * (1 / 𝜇g)
    term3_numerator = f"{c} * {pt}"
    term3_denominator = f"{mu_g}"
    
    # R0f is the product of these three terms. We combine the numerators and denominators.
    final_numerator = f"{term1_numerator} * {term2_numerator} * {term3_numerator}"
    final_denominator = f"{term1_denominator} * {term2_denominator} * {term3_denominator}"
    
    print("The derived expression for R0f is:")
    # We rearrange for clarity, grouping single variables first.
    print(f"R0f = ({b} * {c} * {pg} * {pt} * {tau}) / (({gamma_t} + {mu_t}) * ({tau} + {mu_g}) * {mu_g})")

solve_fire_spread_expression()