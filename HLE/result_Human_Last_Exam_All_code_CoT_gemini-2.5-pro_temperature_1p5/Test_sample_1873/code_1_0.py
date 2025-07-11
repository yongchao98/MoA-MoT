import math
from scipy.stats import norm

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio 'r' based on knowledge graph properties.
    """
    # Given parameters from the problem description
    gamma = 2.1  # Power-law exponent
    alpha = 2.5  # Pareto shape for neighborhood similarity
    epsilon = 0.05 # Marginal completeness tolerance
    confidence = 0.99 # Confidence level

    # Step 1: Calculate the Z-score for the 99% confidence level
    # This corresponds to the 99.5th percentile of a standard normal distribution.
    z_score = norm.ppf(1 - (1 - confidence) / 2)

    # Step 2: Use the most conservative value for completeness p=0.5 to maximize variance.
    p = 0.5
    p_variance = p * (1 - p)

    # Step 3: Define terms for the model r/(1-r) = U_stat / C_graph
    # U_stat represents the squared statistical precision requirement per unit of variance.
    # C_graph represents the structural simplicity of the graph.
    
    # Calculate the statistical uncertainty term
    u_stat = (epsilon / z_score)**2 / p_variance

    # Calculate the graph simplicity term.
    # This term is small for complex graphs (gamma near 2.0).
    c_graph = (gamma - 2) * (alpha - 1)
    
    # Step 4: Calculate the ratio K = U_stat / C_graph
    if c_graph <= 0:
        raise ValueError("Graph complexity term must be positive.")
    k = u_stat / c_graph

    # Step 5: Solve for r from r/(1-r) = K
    r = k / (1 + k)

    # Step 6: Print the calculation step-by-step
    print("Step 1: Determine the Z-score for a confidence level of {:.2f}".format(confidence))
    print("  Z = {:.4f}".format(z_score))
    print("\nStep 2: Define statistical and graph complexity terms")
    print("  Statistical Uncertainty (U_stat) = ((ε/Z)² / (p*(1-p)))")
    print("  U_stat = (({:.2f}/{:.4f})² / {:.2f}) = {:.8f}".format(epsilon, z_score, p_variance, u_stat))
    print("  Graph Simplicity (C_graph) = (γ-2) * (α-1)")
    print("  C_graph = ({:.1f}-2) * ({:.1f}-1) = {:.4f}".format(gamma, alpha, c_graph))
    print("\nStep 3: Calculate the ratio K = U_stat / C_graph")
    print("  K = {:.8f} / {:.4f} = {:.6f}".format(u_stat, c_graph, k))
    print("\nStep 4: Solve for the sampling ratio r from r = K / (1 + K)")
    print("  r = {:.6f} / (1 + {:.6f}) = {:.8f}".format(k, k, r))
    print("\nFinal Answer: The minimum required sampling ratio r, rounded to 4 places, is:")
    print("{:.4f}".format(r))
    
    return r

# Execute the calculation and capture the final value for the specified answer format.
final_r = calculate_sampling_ratio()
# The required output format is <<<answer>>>
# final_answer_str = "<<<{:.4f}>>>".format(final_r)
# print(final_answer_str) # This line would be for a different system. We print the value above.
final_answer_for_submission = round(final_r, 4)
# print(final_answer_for_submission)