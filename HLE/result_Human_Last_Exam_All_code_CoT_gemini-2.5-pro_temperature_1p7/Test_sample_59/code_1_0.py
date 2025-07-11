import scipy.integrate as spi

# According to the Aldous-Hoover theorem for jointly exchangeable graphs,
# the probability of an edge P(y_ij = 1) is determined by a symmetric
# function f(u,v), called a graphon, integrated over the unit square.
# The latent variables u and v are drawn from a uniform distribution U(0,1).

# Since the specific graphon function f(u,v) is not provided, we will use
# a simple, illustrative example to demonstrate the calculation.
# Let's choose the example graphon: f(u, v) = u * v.
def example_graphon(u, v):
  """A simple example of a graphon function f(u,v) = u * v."""
  return u * v

# The integration is over the unit square, so the limits for both u and v
# are from 0 to 1.
lower_bound = 0
upper_bound = 1

# We calculate the double integral of our example graphon to find the probability.
# The dblquad function returns a tuple (result, error_estimate). We are
# interested in the result.
probability, _ = spi.dblquad(
    example_graphon,
    lower_bound,
    upper_bound,
    lambda x: lower_bound, # Inner integral (dv) lower bound
    lambda x: upper_bound  # Inner integral (dv) upper bound
)

# --- Output ---
# First, print the general formula for the probability.
print("The probability P(y_ij = 1) for a jointly exchangeable graph is given by the general formula:")
print("P(y_ij = 1) = ∫ [from u=0 to 1] ∫ [from v=0 to 1] f(u,v) dv du")
print("\nwhere f(u,v) is the graphon defining the graph model.")

# Next, print the specific calculation for our chosen example.
print("\nFor our concrete example using the graphon f(u,v) = u * v, the final equation is:")
# The following print statement outputs each number (0, 1) in the equation.
print(f"P(y_ij = 1) = ∫ (from u={lower_bound} to {upper_bound}) ∫ (from v={lower_bound} to {upper_bound}) (u * v) dv du")

# Finally, print the numerically calculated probability for the example.
print(f"\nThe calculated probability for this specific example is: {probability}")
