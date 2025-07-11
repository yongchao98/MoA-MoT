# The firing level of the rule's antecedent (the "IF" part)
phi_k = 0.7

# The membership value of the rule's consequent (the "THEN" part)
mu_G = 0.9

# The rule activation level is calculated using a t-norm operation, which
# models a fuzzy "AND". The most common t-norm is the minimum operator.
activation_level = min(phi_k, mu_G)

print("Calculating the rule activation level using the minimum t-norm:")
print(f"Antecedent Firing Level Phi_k(x') = {phi_k}")
print(f"Consequent Membership Value μG_k_j(y_j) = {mu_G}")
print("\nThe activation level (w) is the t-norm of these two values:")
print(f"w = min(Phi_k(x'), μG_k_j(y_j))")
print(f"w = min({phi_k}, {mu_G})")
print(f"w = {activation_level}")