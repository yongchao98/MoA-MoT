import scipy.stats
import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.
    
    The plan is as follows:
    1. Define the given parameters: confidence level, tolerance (epsilon),
       power-law exponent (gamma), and Pareto shape (alpha).
    2. Calculate the Z-score for the given confidence level.
    3. Calculate the design effect (deff) which quantifies the increase in variance
       due to the complex graph structure. We model this as the sum of contributions
       from the graph's scale-free nature (gamma) and the neighborhood similarity
       distribution (alpha). The contribution of each is modeled as 1/(exponent - 2).
    4. Formulate and calculate the sampling ratio 'r'. This ratio is modeled as a function
       of the design effect and the statistical requirements (epsilon and Z-score).
       The formula used is r = deff * (epsilon / Z)^2.
    5. Print the final equation with all the intermediate values and the final result.
    """
    
    # 1. Define parameters
    confidence_level = 0.99
    epsilon = 0.05
    gamma = 2.1
    alpha = 2.5
    
    # 2. Calculate Z-score for the two-sided confidence interval
    # alpha_risk is the total probability in the tails (1 - confidence)
    alpha_risk = 1 - confidence_level
    # Z-score corresponds to the point where the cumulative probability is 1 - alpha_risk/2
    z_score = scipy.stats.norm.ppf(1 - alpha_risk / 2)
    
    # 3. Calculate the design effect (deff)
    # The design effect is modeled as the sum of complexity factors from gamma and alpha.
    deff = (1 / (gamma - 2)) + (1 / (alpha - 2))
    
    # 4. Calculate the minimum sampling ratio (r)
    # The ratio r is modeled as deff * (epsilon / z_score)^2
    r = deff * (epsilon / z_score)**2
    
    # 5. Print the final calculation and result
    print(f"Z-score for {confidence_level*100}% confidence: {z_score:.4f}")
    print(f"Marginal completeness tolerance (ε): {epsilon}")
    print(f"Power-law exponent (γ): {gamma}")
    print(f"Pareto shape (α): {alpha}")
    print("\nCalculating the design effect (deff):")
    print(f"deff = 1/(γ - 2) + 1/(α - 2)")
    print(f"deff = 1/({gamma} - 2) + 1/({alpha} - 2)")
    print(f"deff = 1/{gamma - 2:.1f} + 1/{alpha - 2:.1f}")
    print(f"deff = {1/(gamma - 2):.4f} + {1/(alpha - 2):.4f} = {deff:.4f}")
    
    print("\nCalculating the minimum sampling ratio (r):")
    print(f"r = deff * (ε / Z)²")
    print(f"r = {deff:.4f} * ({epsilon} / {z_score:.4f})²")
    print(f"r = {deff:.4f} * ({epsilon / z_score:.4f})²")
    print(f"r = {deff:.4f} * {(epsilon / z_score)**2:.8f}")
    print(f"r = {r:.4f}")

calculate_sampling_ratio()
<<<0.0045>>>