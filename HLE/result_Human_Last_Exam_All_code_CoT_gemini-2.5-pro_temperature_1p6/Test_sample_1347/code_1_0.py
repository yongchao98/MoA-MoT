import math

def solve_fire_spread():
    """
    This script derives the expression for R0f and demonstrates its calculation.
    """
    
    # Part 1: Explain the derivation of the R0f expression
    explanation = """
The expression for R0f, the number of additional trees that catch fire from a single burning tree, is derived by analyzing the two-step infection process: Tree -> Grass -> Tree.

Step 1: Grass patches ignited by a single tree.
A single burning tree is infectious for an average period of 1 / (𝛾t + 𝜇t), where 𝛾t is the extinguishing rate and 𝜇t is the tree's natural death rate.
During this time, it successfully ignites grass at a rate of (b * pg).
So, the total number of grass patches ignited by one tree is:
Component A = (b * pg) / (𝛾t + 𝜇t)

Step 2: Probability of ignited grass becoming infectious.
Ignited grass enters a latent period. It becomes infectious at a rate 𝜏, but can also die at a rate 𝜇g. The probability of it surviving this period to become infectious is the ratio of the success rate to the total rate of leaving the latent state:
Component B = 𝜏 / (𝜏 + 𝜇g)

Step 3: Trees ignited by a single infectious grass patch.
An infectious grass patch spreads fire for an average period of 1 / 𝜇g (its lifetime until natural death).
During this time, it successfully ignites trees at a rate of (c * pt).
So, the total number of trees ignited by one infectious grass patch is:
Component C = (c * pt) / 𝜇g

Final Expression for R0f:
R0f is the product of these three components: A * B * C
R0f = [ (b * pg) / (𝛾t + 𝜇t) ] * [ 𝜏 / (𝜏 + 𝜇g) ] * [ (c * pt) / 𝜇g ]

This simplifies to:
R0f = (b * pg * c * pt * 𝜏) / ( (𝛾t + 𝜇t) * (𝜏 + 𝜇g) * 𝜇g )
"""
    print(explanation)

    # Part 2: Define a Python function for the formula
    def calculate_r0f(b, pg, c, pt, gamma_t, mu_t, tau, mu_g):
        """
        Calculates R0f based on the derived formula.
        """
        numerator = b * pg * c * pt * tau
        denominator = (gamma_t + mu_t) * (tau + mu_g) * mu_g
        if denominator == 0:
            return float('inf')
        return numerator / denominator

    # Part 3: Assign and print example values
    # These are hypothetical values for a clear demonstration
    b = 5.0      # number of grasses each burning tree can ignite per day
    pg = 0.1     # probability that contact (tree->grass) ignites grass
    c = 10.0     # number of trees each smoldering area can ignite per day
    pt = 0.05    # probability that contact (grass->tree) ignites tree
    gamma_t = 0.1 # fire extinguishing rate for burning trees (per day)
    mu_t = 0.01  # rate at which trees naturally die (per day)
    tau = 2.0    # rate ignited grass starts burning intensely (1/𝜏 is the avg time)
    mu_g = 0.2   # per-day rate at which grass naturally dies

    print("\n--- Example Calculation ---")
    print("Using the following example parameter values:")
    print(f"b = {b}, pg = {pg}, c = {c}, pt = {pt}, 𝛾t = {gamma_t}, 𝜇t = {mu_t}, 𝜏 = {tau}, 𝜇g = {mu_g}\n")
    
    # Part 4: Calculate the result with the example values
    r0f_value = calculate_r0f(b, pg, c, pt, gamma_t, mu_t, tau, mu_g)

    # Part 5: Print the final equation with each number substituted in
    print("The final equation with these numbers is:")
    
    numerator_str = f"({b} * {pg} * {c} * {pt} * {tau})"
    denominator_str = f"(({gamma_t} + {mu_t}) * ({tau} + {mu_g}) * {mu_g})"
    print(f"R0f = {numerator_str} / {denominator_str}")

    # Calculate and show intermediate values for clarity
    numerator_val = b * pg * c * pt * tau
    term1_denom_val = gamma_t + mu_t
    term2_denom_val = tau + mu_g
    term3_denom_val = mu_g
    denominator_val = term1_denom_val * term2_denom_val * term3_denom_val
    
    print(f"R0f = {numerator_val} / ({term1_denom_val:.2f} * {term2_denom_val:.1f} * {term3_denom_val})")
    print(f"R0f = {numerator_val} / {denominator_val:.4f}")
    print(f"R0f = {r0f_value:.4f}")

solve_fire_spread()