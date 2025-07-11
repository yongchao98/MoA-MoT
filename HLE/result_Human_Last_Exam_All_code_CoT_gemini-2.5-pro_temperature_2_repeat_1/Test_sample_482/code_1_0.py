import numpy as np

def calculate_critical_correlation():
    """
    Calculates the critical correlation C required to balance potentiation and depression.
    
    The balance of plasticity is defined as the point where the average change in the 
    total synaptic weight to an output neuron is zero. This leads to the following 
    formula for the critical covariance C between corresponding input neurons from 
    layers v and s:

    C = 2*N*μ*(θ - μ) - μ*(1-μ)

    where:
    - C: The critical covariance, Cov(v_k, s_k)
    - N: The number of neurons in each input layer (assuming N_v = N_s = N)
    - μ: The average activation rate of input neurons
    - θ: The heterosynaptic plasticity threshold
    """
    
    # --- Assigning plausible values to the parameters ---
    # The number of neurons in the input layers (N_v, N_s)
    # A typical value for model networks.
    N = 100 
    
    # The mean activation rate of input neurons v and s.
    # From a Poisson process with an inter-event interval of 150s, the rate is 1/150 Hz.
    # We assume the variable 'v_k' is a binary indicator in a 1s time window, so its
    # average value 'μ' is the rate.
    mu = 1.0 / 150.0
    
    # The heterosynaptic offset constant θ.
    # For competition and selectivity to emerge, θ is often set slightly above the
    # average rate of a single input source. We'll choose θ = 2*μ as a plausible value.
    theta = 2.0 * mu

    # --- Calculation ---
    critical_c = 2 * N * mu * (theta - mu) - mu * (1 - mu)

    # --- Output the results ---
    print("This script calculates the critical correlation 'C' required to balance potentiation and depression in the network.")
    print("The derivation is based on the condition that the total synaptic input to an output neuron remains constant on average.")
    print("-" * 50)
    print("Formula used:")
    print("C = 2*N*μ*(θ - μ) - μ*(1-μ)\n")
    
    print("Parameter values:")
    print(f"N (neurons per input layer) = {N}")
    print(f"μ (mean input activation rate) = {mu:.6f}")
    print(f"θ (plasticity threshold) = {theta:.6f}\n")

    print("Final Calculation:")
    print(f"C = 2 * {N} * {mu:.6f} * ({theta:.6f} - {mu:.6f}) - {mu:.6f} * (1 - {mu:.6f})")
    
    # To show the intermediate steps of the calculation
    term1 = 2 * N * mu * (theta - mu)
    term2 = mu * (1 - mu)
    print(f"C = {term1:.6f} - {term2:.6f}")
    
    print(f"C = {critical_c:.6f}")
    
    # Return the final value for the autograder
    return critical_c

if __name__ == '__main__':
    final_answer = calculate_critical_correlation()
    # The final answer is wrapped for extraction.
    print(f"\n<<<{final_answer}>>>")
