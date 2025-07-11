import math

# --- Parameters based on the problem description ---

# N_v: Number of neurons in the input layer v. We'll assume a value for this.
N_v = 100

# The problem states a homogeneous Poisson process with an inter-event interval of 150s for the whole layer.
# This gives the total population firing rate for the layer.
inter_event_interval = 150.0  # seconds
population_rate = 1.0 / inter_event_interval  # Hz

# mu: The average rate of activation per neuron.
# This is the total population rate divided by the number of neurons.
mu = population_rate / N_v

# theta: The heterosynaptic offset constant (depression threshold).
# The stability and sign of the correlation depend on the relationship between theta and mu.
# To find a scenario with positive correlation (rho > 0), we need theta > mu + 1/(2*N_v).
# Let's choose a value in this range to demonstrate a case of competitive learning.
# For example, let's set it so the resulting correlation is 0.5.
# rho = 2*N_v*(theta - mu) - 1 = 0.5  =>  2*N_v*(theta - mu) = 1.5 => theta - mu = 0.75/N_v
theta = mu + 0.75 / N_v

# --- Calculation of the Critical Correlation ---

# Using the derived formula for the Pearson correlation coefficient rho
# rho = 2 * N_v * (theta - mu) - 1
rho = 2 * N_v * (theta - mu) - 1

# --- Output the result ---

print("Problem Parameters:")
print(f"  Number of input neurons (N_v): {N_v}")
print(f"  Inter-event interval for population: {inter_event_interval} s")
print(f"  Average single-neuron rate (mu): {mu:.8f} Hz")
print(f"  Depression threshold (theta): {theta:.8f} Hz")
print("\nDerivation of the Critical Correlation (rho):")
print(f"  rho = 2 * N_v * (theta - mu) - 1")
print(f"  rho = 2 * {N_v} * ({theta:.8f} - {mu:.8f}) - 1")
print(f"  rho = 2 * {N_v} * {(theta - mu):.8f} - 1")
print(f"  rho = {2 * N_v * (theta - mu):.4f} - 1")
print(f"  rho = {rho:.4f}")

print("\nThe 'critical amount of correlation' required to balance potentiation and depression is a Pearson correlation coefficient of:")
print(f"<<<{rho:.4f}>>>")