import math

def calculate_critical_correlation():
    """
    Calculates the critical correlation C required to balance potentiation and depression.
    
    The critical correlation is derived from the stability condition of the weight dynamics,
    leading to the formula: C = mu * (1 - mu), where C is the covariance between
    corresponding input neurons (v_k, s_k) and mu is their average firing probability.
    """
    
    # Given parameters
    # Inter-event interval for the Poisson process
    inter_event_interval = 150.0  # seconds

    # Assumed simulation time step (dt). The problem states "sufficiently small time steps".
    # 1ms is a standard choice in computational neuroscience simulations.
    dt = 0.001  # seconds
    
    # --- Calculations ---
    # 1. Calculate the rate (lambda) of the Poisson process
    rate_lambda = 1.0 / inter_event_interval
    
    # 2. Calculate the firing probability (mu) in a single time step dt
    # For a small dt, P(event) in a Poisson process is approximately lambda * dt
    mu = rate_lambda * dt
    
    # 3. Calculate the critical covariance (C)
    critical_covariance = mu * (1 - mu)
    
    # --- Output Results ---
    print("Derivation Summary:")
    print("The condition to balance potentiation and depression while allowing for synaptic selectivity is:")
    print("C = μ * (1 - μ)")
    print("Where C is the covariance Cov(v_k, s_k) and μ is the firing probability per time step.\n")

    print("Calculation with Given and Assumed Values:")
    print(f"  Inter-event interval T = {inter_event_interval} s")
    print(f"  Assumed time step dt = {dt} s")
    print("-" * 20)
    print(f"  Calculated firing rate λ (1/T) = {rate_lambda:.6f} Hz")
    print(f"  Calculated firing probability μ (λ * dt) = {mu:.6e}")
    print(f"  Calculated critical covariance C (μ * (1 - μ)) = {critical_covariance:.6e}\n")

    print("Final Equation with Numbers:")
    # Using f-string formatting to show the full numerical equation
    print(f"{critical_covariance:.6e} = {mu:.6e} * (1 - {mu:.6e})")
    print(f"{critical_covariance:.6e} = {mu:.6e} * {1-mu:.6f}")

    # Returning the final value as requested by the output format
    return critical_covariance

# Run the calculation and store the result
final_answer = calculate_critical_correlation()

# The final answer as a raw value
# <<<final_answer>>>