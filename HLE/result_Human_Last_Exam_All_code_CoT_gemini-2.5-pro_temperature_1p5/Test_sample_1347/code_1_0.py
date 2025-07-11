import sys
from platform import python_version

# Check for a specific Python version or other environment settings if needed
# For this script, we'll just print the expression, so no special checks are required.

def solve_r0f_expression():
    """
    This function constructs and prints the symbolic expression for R0f
    based on the provided parameters of the fire spread model.
    """
    
    # Define the symbols as strings for clarity in the final expression
    b = 'b'
    pg = 'pg'
    tau = '𝜏'
    c = 'c'
    pt = 'pt'
    gamma_t = '𝛾t'
    mu_t = '𝜇t'
    mu_g = '𝜇g'

    # Part 1: Number of grass areas ignited by one tree over its lifetime
    # Lifetime of a burning tree = 1 / (𝛾t + 𝜇t)
    # Rate of grass ignition = b * pg
    grass_ignited_per_tree = f"({b} * {pg}) / ({gamma_t} + {mu_t})"

    # Part 2: Probability of ignited grass becoming intensely burning
    # Rate of becoming intense = 𝜏
    # Rate of removal before becoming intense = 𝜇g
    prob_grass_matures = f"{tau} / ({tau} + {mu_g})"

    # Part 3: Number of trees ignited by one intensely burning grass area over its lifetime
    # Lifetime of burning grass = 1 / 𝜇g
    # Rate of tree ignition = c * pt
    trees_ignited_per_grass = f"({c} * {pt}) / {mu_g}"

    # R0f is the product of these three parts
    # R0f = (Part 1) * (Part 2) * (Part 3)
    
    # Numerator of the final expression
    numerator_expression = f"{b} * {pg} * {tau} * {c} * {pt}"

    # Denominator of the final expression
    denominator_expression = f"({gamma_t} + {mu_t}) * ({tau} + {mu_g}) * {mu_g}"

    print("The expression for R0f is derived by multiplying the contributions from each stage of the fire's spread from a tree to other trees via grass.")
    
    print("\n1. Number of grass areas ignited per tree:", grass_ignited_per_tree)
    print("2. Probability of grass becoming intensely burning:", prob_grass_matures)
    print("3. Number of new trees ignited per intensely burning grass area:", trees_ignited_per_grass)
    
    print("\nCombining these terms gives the final expression for R0f:")
    print(f"R0f = ({numerator_expression}) / ({denominator_expression})")

# Execute the function to print the derivation and result
solve_r0f_expression()
