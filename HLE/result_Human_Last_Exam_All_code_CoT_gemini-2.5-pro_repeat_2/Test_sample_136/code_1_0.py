def solve_torus_probability():
    """
    Calculates the limit of the conditional probability based on local escape probabilities.
    """
    # At t=1, the walk moves from 0 to one of its 4 neighbors with probability 1/4 each.
    # The two common neighbors of 0 and x_0 are c1, c2.
    # The other two neighbors of 0 are v1, v2.
    
    # --- Numerator Calculation ---
    # Probability of avoiding {0, x_0} at the next step, averaged over the first step.
    # If the first step is to a common neighbor (c1 or c2), its neighbors include both 0 and x_0.
    # The degree of any vertex is 4.
    prob_hit_S_from_common_neighbor = 2 / 4
    prob_escape_S_from_common_neighbor = 1 - prob_hit_S_from_common_neighbor

    # If the first step is to a non-common neighbor (v1 or v2), its neighbors include 0 but not x_0.
    prob_hit_S_from_non_common_neighbor = 1 / 4
    prob_escape_S_from_non_common_neighbor = 1 - prob_hit_S_from_non_common_neighbor

    # Average escape probability for the numerator (2 common, 2 non-common neighbors)
    avg_escape_numerator = (1/4) * (prob_escape_S_from_common_neighbor * 2 + prob_escape_S_from_non_common_neighbor * 2)

    # --- Denominator Calculation ---
    # Probability of avoiding {0} at the next step, averaged over the first step.
    # For any neighbor y of 0, the probability of returning to 0 is 1/4.
    prob_hit_S_prime = 1 / 4
    prob_escape_S_prime = 1 - prob_hit_S_prime
    
    # Average escape probability for the denominator is the same from all 4 neighbors.
    avg_escape_denominator = (1/4) * (prob_escape_S_prime * 4)

    # --- Final Result ---
    limit = avg_escape_numerator / avg_escape_denominator
    
    # Print the step-by-step calculation
    num_val_1 = prob_escape_S_from_common_neighbor
    num_val_2 = prob_escape_S_from_non_common_neighbor
    print(f"Numerator's average local escape probability: (1/4) * ({num_val_1} + {num_val_1} + {num_val_2} + {num_val_2}) = {avg_escape_numerator}")

    den_val = prob_escape_S_prime
    print(f"Denominator's average local escape probability: (1/4) * ({den_val} + {den_val} + {den_val} + {den_val}) = {avg_escape_denominator}")
    
    print(f"The limit is the ratio: ({avg_escape_numerator}) / ({avg_escape_denominator}) = {limit}")


solve_torus_probability()
<<<5/6>>>