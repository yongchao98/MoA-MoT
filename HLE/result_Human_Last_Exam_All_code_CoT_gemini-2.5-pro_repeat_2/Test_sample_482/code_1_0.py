import math

# --- Parameters ---
# Number of neurons in each input population (v and s)
N = 100
# Average firing rate of input neurons (in Hz)
# Derived from the given inter-event interval of 150s (mu = 1/150)
mu = 1.0 / 150.0
# Heterosynaptic offset constant (in Hz)
# We assume a value slightly larger than mu to demonstrate.
theta = 0.01
# Variance of the input firing rates.
# For a Poisson process, variance of the spike count in time T is mu*T.
# The rate variance is (mu*T)/T^2 = mu/T. We assume T=1s for this example.
sigma_sq = mu 

# --- Calculation ---
# The formula for the critical correlation C is derived from the balance
# between potentiation and depression across the network.
# C = 2 * N * mu * (theta - mu) - sigma^2
critical_correlation = 2 * N * mu * (theta - mu) - sigma_sq

# --- Output ---
print("Calculating the critical amount of correlation (C)")
print("-" * 50)
print(f"Parameters used:")
print(f"  N (neurons per population) = {N}")
print(f"  mu (average input rate)    = {mu:.6f} Hz")
print(f"  theta (depression offset)  = {theta:.6f} Hz")
print(f"  sigma^2 (input variance)   = {sigma_sq:.6f} Hz^2")
print("-" * 50)
print("Formula: C = 2 * N * mu * (theta - mu) - sigma^2")
print(f"Calculation: C = 2 * {N} * {mu:.6f} * ({theta:.6f} - {mu:.6f}) - {sigma_sq:.6f}")
print("-" * 50)
print(f"The critical correlation C is: {critical_correlation:.6f}")

# The final answer in the required format
# print(f"<<<{critical_correlation}>>>")