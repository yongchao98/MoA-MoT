import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # Step 1: Define all parameters based on the problem description.
    # We use milliseconds (ms) for time and millivolts (mV) for potential.
    tau_m = 20.0      # Membrane time constant (ms)
    V_reset = 10.0    # Voltage reset (mV)
    V_th = 20.0       # Voltage threshold (mV)
    tau_ref = 2.0     # Refractory period (ms)
    g = 4.0           # Relative inhibition to excitation factor
    K_E = 1000.0      # Number of excitatory connections
    V_ext = 30.0      # External input (mV)

    # Step 2: Calculate the number of inhibitory connections.
    # The problem states "a 0.25 proportion of that", where "that" refers to
    # the number of excitatory connections (K_E).
    K_I = 0.25 * K_E

    # Step 3: Calculate the mean membrane potential (μ).
    # The mean potential of a neuron in the network is determined by synaptic
    # inputs and external drive. In this specific configuration, the network
    # is balanced. Let's verify:
    balance_term = K_E - g * K_I
    
    # Because the balance term is 0, the mean synaptic input is zero,
    # and the mean membrane potential is simply equal to the external input.
    mu = V_ext

    # Step 4: Calculate the neuron's firing rate using the LIF formula.
    firing_rate_hz = 0
    if mu > V_th:
        # The time between spikes is the sum of the refractory period and the time
        # it takes for the membrane potential to rise from V_reset to V_th.
        # Firing Rate (in Hz) = 1000 / (time between spikes in ms)
        time_between_spikes_ms = tau_ref + tau_m * math.log((mu - V_reset) / (mu - V_th))
        firing_rate_hz = 1000.0 / time_between_spikes_ms

    # Step 5: Print the explanation and the final result.
    print("--- Firing Rate Calculation Steps ---")
    print(f"\n1. The network is in a balanced state because K_E - g * K_I = {K_E} - {g} * {K_I} = {balance_term}.")
    print(f"   Therefore, the mean membrane potential (μ) is determined solely by the external input.")
    print(f"   μ = {mu:.1f} mV")

    print(f"\n2. Since μ ({mu:.1f} mV) is greater than V_th ({V_th:.1f} mV), the neuron will fire.")
    print("   The firing rate (ν) is calculated using the formula:")
    print("   ν(Hz) = 1000 / (τ_ref + τ * ln((μ - V_reset) / (μ - V_th)))")

    # Print the final equation with all the numerical values plugged in.
    print("\nFinal Equation:")
    final_equation = f"ν(Hz) = 1000 / ({tau_ref} + {tau_m} * math.log(({mu} - {V_reset}) / ({mu} - {V_th})))"
    print(final_equation)

    # Print the final numerical result.
    result_int = int(round(firing_rate_hz))
    print(f"\nCalculated firing rate: {firing_rate_hz:.2f} Hz")
    print(f"\nThe firing rate of a typical neuron is {result_int} Hz.")
    
    return result_int

# Execute the calculation and store the final integer answer.
final_answer = calculate_firing_rate()
print(f"<<<{final_answer}>>>")