import numpy as np

def calculate_critical_correlation():
    """
    This function calculates the critical correlation 'c' required to balance
    potentiation and depression in the described neural network.

    The balance is achieved when the total synaptic weight onto an output neuron is stable.
    This leads to the following derived formula for the critical correlation c:
    c = (N_v + N_s) * (θ - μ) - 1
    
    where:
    - N_v, N_s: Number of neurons in input layers v and s.
    - θ: Heterosynaptic depression constant.
    - μ: Average rate of activation for input neurons.
    - c: The correlation coefficient between corresponding v and s neurons.
    """

    # Since no specific parameters were given, we use some plausible example values.
    N_v = 50  # Number of neurons in input layer v
    N_s = 50  # Number of neurons in input layer s
    mu = 0.05 # Average firing rate of input neurons (e.g., in Hz)
    theta = 0.06 # Plasticity threshold (e.g., in Hz)

    # Calculate the critical correlation using the derived formula.
    # c = (N_v + N_s) * (theta - mu) - 1
    term1 = N_v + N_s
    term2 = theta - mu
    critical_c = term1 * term2 - 1
    
    # Print the explanation and the result.
    print("The final equation for the critical correlation 'c' is:")
    print("c = (N_v + N_s) * (θ - μ) - 1\n")
    print("Plugging in the example values:")
    print(f"N_v = {N_v}")
    print(f"N_s = {N_s}")
    print(f"θ   = {theta}")
    print(f"μ   = {mu}\n")
    print("The equation becomes:")
    # Using 'g' format to avoid unnecessary trailing zeros for integers
    print(f"c = ({N_v:g} + {N_s:g}) * ({theta:g} - {mu:g}) - 1")
    print(f"c = {term1:g} * {term2:g} - 1")
    print(f"c = {term1 * term2:g} - 1\n")
    print(f"Calculated critical correlation c = {critical_c:.4f}")

if __name__ == '__main__':
    calculate_critical_correlation()
