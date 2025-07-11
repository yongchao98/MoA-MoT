import math

def calculate_critical_correlation(Nv, Ns, mu, theta, tau_r):
    """
    Calculates the critical amount of correlation 'C' required to balance potentiation and depression.

    The balance is achieved when the expected weight change is zero. This leads to the formula:
    C = mu * ((Nv + Ns) * (theta - mu) - (1 / tau_r))

    Args:
        Nv (int): Number of neurons in the input layer v.
        Ns (int): Number of neurons in the input layer s.
        mu (float): Average rate of activation for v and s (in Hz).
        theta (float): Heterosynaptic offset constant (in Hz).
        tau_r (float): Leaky integrator time constant for r (in seconds).
    
    Returns:
        float: The critical correlation C (Cov(v,s), in Hz^2).
    """

    print("Calculating the critical correlation C...")
    print(f"The formula derived is: C = mu * ((Nv + Ns) * (theta - mu) - 1/tau_r)")
    print("-" * 30)
    print("Parameter values:")
    print(f"  Number of neurons Nv: {Nv}")
    print(f"  Number of neurons Ns: {Ns}")
    print(f"  Average input rate mu: {mu:.4f} Hz")
    print(f"  Heterosynaptic offset theta: {theta:.4f} Hz")
    print(f"  Time constant tau_r: {tau_r:.4f} s")
    print("-" * 30)

    # To be physically plausible, the covariance C must satisfy |C| <= Var(v).
    # Assuming Var(v) = mu / tau_r, we must choose theta such that:
    # mu <= theta <= mu + 2 / (tau_r * (Nv + Ns))
    var_v = mu / tau_r
    max_theta = mu + 2 / (tau_r * (Nv + Ns))
    print(f"For these parameters, the variance Var(v) is assumed to be mu/tau_r = {var_v:.4f} Hz^2.")
    print(f"For a valid correlation, theta must be in [{mu:.4f}, {max_theta:.4f}]. The chosen theta = {theta:.4f} is in this range.")
    print("-" * 30)
    
    # Calculate each part of the equation
    term1 = Nv + Ns
    term2 = theta - mu
    term3 = 1 / tau_r
    
    print("Calculation steps:")
    print(f"  (Nv + Ns) = {term1}")
    print(f"  (theta - mu) = {theta:.4f} - {mu:.4f} = {term2:.4f}")
    print(f"  1 / tau_r = 1 / {tau_r:.4f} = {term3:.4f}")
    
    critical_c = mu * (term1 * term2 - term3)

    print("\nFinal calculation:")
    print(f"  C = {mu:.4f} * (({term1}) * ({term2:.4f}) - ({term3:.4f}))")
    print(f"  C = {mu:.4f} * ({term1 * term2:.4f} - {term3:.4f})")
    print(f"  C = {mu:.4f} * ({term1 * term2 - term3:.4f})")
    print(f"  C = {critical_c:.4f} Hz^2")

    return critical_c

# --- Example Usage ---
# Problem states inter-event interval of 150s for the Poisson process.
# This corresponds to a mean firing rate mu.
mu_rate = 1.0 / 150.0  # in Hz

# Let's define the other parameters with reasonable values.
# Number of neurons in each input layer
N_v_neurons = 100
N_s_neurons = 100

# Leaky integrator time constant (e.g., 20ms)
tau_r_constant = 0.02  # in seconds

# Heterosynaptic offset constant, chosen to be in the valid range.
# For these values, theta should be between ~0.0067 and 0.2567 Hz.
# Let's pick a value in the middle of this range.
theta_constant = 0.15

# Calculate and print the result
critical_correlation = calculate_critical_correlation(
    Nv=N_v_neurons,
    Ns=N_s_neurons,
    mu=mu_rate,
    theta=theta_constant,
    tau_r=tau_r_constant
)
print(f"\n<<<The critical amount of correlation is {critical_correlation:.4f} Hz^2>>>")