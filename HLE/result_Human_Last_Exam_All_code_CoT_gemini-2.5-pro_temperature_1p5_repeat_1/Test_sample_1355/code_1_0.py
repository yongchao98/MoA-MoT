import math

# Explain the theoretical background and the formula used.
print(
    "For an ensemble of disordered Majorana wires at the critical point, the statistics of the\n"
    "dimensionless conductance 'g' are universal in the limit of a large wire length. The n-th\n"
    "statistical moment of the conductance is given by a specific formula from random-matrix theory."
)
print("\nThe formula for the n-th moment, <g^n>, is:")
print("  <g^n> = Γ(n + 1/2) / (Γ(n + 1) * Γ(1/2))")
print("where Γ(z) is the Euler Gamma function.\n")

# Step 1: Calculate the average conductance <g> (the 1st moment, n=1).
n_avg = 1
avg_g = math.gamma(n_avg + 0.5) / (math.gamma(n_avg + 1) * math.gamma(0.5))
print(f"Step 1: Calculate the average value <g> using n = {n_avg}.")
print(f"  <g> = {avg_g}\n")

# Step 2: Calculate the fourth statistical moment <g^4> (n=4).
n_fourth = 4
fourth_moment_g = math.gamma(n_fourth + 0.5) / (math.gamma(n_fourth + 1) * math.gamma(0.5))
print(f"Step 2: Calculate the fourth moment <g^4> using n = {n_fourth}.")
print(f"  <g^4> = {fourth_moment_g}\n")

# Step 3: Calculate the final ratio.
ratio = fourth_moment_g / avg_g
print("Step 3: Compute the ratio between the fourth moment and the average value.")

# Display the final result in the form of an equation.
print("\nThe final equation is:")
print(f"Ratio = <g^4> / <g> = {fourth_moment_g} / {avg_g} = {ratio}")
