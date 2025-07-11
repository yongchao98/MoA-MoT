import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # Step 1: Define the network and neuron parameters.
    # All time units are in ms, and voltage units are in mV.
    tau_m = 20.0      # Membrane time constant (ms)
    V_reset = 10.0    # Voltage reset (mV)
    V_th = 20.0       # Voltage threshold (mV)
    tau_ref = 2.0     # Refractory period (ms)
    V_ext = 30.0      # External input (mV)
    K_E = 1000.0      # Number of excitatory connections
    p_I = 0.25        # Proportion of inhibitory to excitatory connections
    g = 4.0           # Relative inhibition to excitation strength

    # Step 2: Analyze the average synaptic input.
    # The network is in a balanced state where mean excitatory and inhibitory inputs cancel.
    # K_I = p_I * K_E = 0.25 * 1000 = 250
    # The net effect of synaptic input is proportional to (K_E - g * K_I).
    # (1000 - 4 * 250) = 1000 - 1000 = 0.
    # Therefore, the average synaptic input is zero.
    
    # Step 3: Determine the total mean input (mu).
    # Since the average synaptic input is zero, the total mean input is the external input.
    mu = V_ext

    print(f"The total mean input to the neuron (μ) is equal to the external input: {mu} mV.")
    print(f"This is because the network is in a balanced state where excitatory and inhibitory inputs cancel out on average.")
    print("-" * 50)
    
    # Check if the neuron is in a supra-threshold regime.
    if mu <= V_th:
        print("The mean input is not above the threshold, so the firing rate would be 0 Hz in this simplified model.")
        return

    # Step 4: Use the LIF neuron firing rate formula.
    # The formula for the inter-spike interval (T_fire) is:
    # T_fire = tau_ref + tau_m * ln((mu - V_reset) / (mu - V_th))
    
    # Calculate terms for the equation.
    numerator = mu - V_reset
    denominator = mu - V_th
    
    # Step 5: Calculate the result and print the final equation with numbers.
    print("Calculating the firing rate using the formula for a leaky integrate-and-fire neuron.")
    print("The final equation for the inter-spike interval (T) with the given values is:")
    
    # Print the equation with all numbers substituted.
    final_equation = f"T = {tau_ref} + {tau_m} * ln(({mu} - {V_reset}) / ({mu} - {V_th}))"
    print(final_equation)
    
    simplified_equation = f"T = {tau_ref} + {tau_m} * ln({numerator} / {denominator})"
    print(simplified_equation)
    
    # Perform the calculation.
    T_fire_ms = tau_ref + tau_m * math.log(numerator / denominator)
    
    # Firing rate in Hz is 1000 / T_fire (in ms)
    firing_rate_hz = 1000 / T_fire_ms
    
    # Round to the nearest integer as requested.
    final_rate = int(round(firing_rate_hz))
    
    print("-" * 50)
    print(f"The calculated inter-spike interval is {T_fire_ms:.2f} ms.")
    print(f"The corresponding firing rate is {firing_rate_hz:.2f} Hz.")
    print(f"\nThe firing rate of a typical neuron, rounded to the nearest integer, is: {final_rate}")


if __name__ == "__main__":
    calculate_firing_rate()
    print("\n<<<63>>>")