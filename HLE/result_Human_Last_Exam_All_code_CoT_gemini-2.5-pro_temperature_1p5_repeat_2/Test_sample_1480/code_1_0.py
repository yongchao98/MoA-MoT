# Define the given membership values from the problem statement.
phi_val = 0.7
mu_val = 0.9

# The rule activation level is calculated by applying a t-norm to the antecedent
# membership values. We will use the minimum t-norm, which is the most common choice.
# The calculation is the minimum of the two provided values.
activation_level = min(phi_val, mu_val)

# Print the final result, showing the full equation as requested.
print(f"The problem is to find the rule activation level using a t-norm on {phi_val} and {mu_val}.")
print(f"Using the minimum t-norm, the calculation is:")
print(f"min({phi_val}, {mu_val}) = {activation_level}")