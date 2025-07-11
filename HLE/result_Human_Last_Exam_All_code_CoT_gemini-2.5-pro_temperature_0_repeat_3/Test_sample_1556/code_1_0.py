# Define the measured probability of obtaining the |0> state at the output.
p0_out = 0.36

# Based on the analysis, the quantum circuit B effectively swaps the probabilities
# of the input state's components.
# Therefore, the probability of measuring |0> at the output is equal to the
# probability of the input state being in the |1> state.
# P(0)_out = P(1)_in = |β|²
beta_squared = p0_out

# The total probability of the input state must be 1.
# |α|² + |β|² = 1
# We can now solve for |α|², which is the probability of the input state being in the |0> state.
alpha_squared = 1 - beta_squared

# Print the equation and the final result.
print(f"The value of |α|² is calculated from the normalization condition: |α|² + |β|² = 1")
print(f"Given P(0)_out = |β|² = {p0_out}")
print(f"Therefore, |α|² = 1 - {beta_squared} = {alpha_squared}")
