# Define the given values for the fuzzy rule inference
# Firing strength of the rule's antecedent (IF part)
antecedent_firing_strength = 0.7
# Membership grade in the consequent set (THEN part)
consequent_membership_grade = 0.9

# A t-norm is used to combine the antecedent's firing strength with the consequent's membership function.
# This process determines the membership grade in the final output set for this rule.
# The most widely used t-norm is the minimum operator.
result = min(antecedent_firing_strength, consequent_membership_grade)

# Display the inputs, the operation, and the final result
print("Calculating the resulting membership grade using a t-norm (minimum operator).")
print(f"Antecedent Firing Strength, Phi_k(x'): {antecedent_firing_strength}")
print(f"Consequent Membership Grade, μG_k_j(y_j): {consequent_membership_grade}")
print("\nEquation:")
print(f"Result = min({antecedent_firing_strength}, {consequent_membership_grade})")
print(f"\nFinal Result: {result}")