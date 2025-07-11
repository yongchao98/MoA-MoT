import math

def calculate_neuron_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network based on given parameters.
    """
    # 1. Define model parameters from the problem description.
    # Time parameters are converted from milliseconds to seconds.
    tau_m = 20.0 / 1000.0   # Membrane time constant in seconds
    tau_ref = 2.0 / 1000.0  # Refractory period in seconds
    V_reset = 10.0          # Reset potential in mV
    V_th = 20.0             # Firing threshold in mV
    mu_ext = 30.0           # External input drive in mV
    g = 4.0                 # Relative strength of inhibition to excitation
    K_E = 1000.0            # Number of excitatory connections
    
    # The number of inhibitory connections is 0.25 times the number of excitatory ones.
    K_I = 0.25 * K_E

    # 2. Determine the total input potential (mu).
    # The net effect of recurrent connections is proportional to (K_E - g * K_I).
    # If this term is zero, the network is balanced and recurrent inputs cancel out.
    balance_factor = K_E - g * K_I
    
    # In a balanced network, the total input potential 'mu' is equal to the external input.
    if balance_factor == 0:
        mu = mu_ext
    else:
        # For this problem, the network is designed to be balanced.
        # This part handles the general case but won't be triggered.
        # Solving for a self-consistent firing rate 'r' would be needed otherwise.
        mu = mu_ext 

    print("Step 1: Determine the total input potential 'mu'.")
    print(f"The analysis shows the network is balanced because the net recurrent input factor (K_E - g * K_I) is {int(K_E)} - {g} * {int(K_I)} = {balance_factor}.")
    print("Therefore, the total input potential 'mu' is simply the external input.")
    print(f"mu = {mu} mV\n")

    # 3. Calculate the firing rate using the LIF formula.
    # The formula is valid because mu (30 mV) > V_th (20 mV).
    if mu > V_th:
        # Components for the formula
        log_arg_numerator = mu - V_reset
        log_arg_denominator = mu - V_th
        log_argument = log_arg_numerator / log_arg_denominator
        
        # Calculate firing rate in Hz
        firing_rate_hz = 1 / (tau_ref + tau_m * math.log(log_argument))
        
        # Round to the nearest integer for the final answer
        final_answer = int(round(firing_rate_hz))
        
        print("Step 2: Calculate the firing rate 'r' using the LIF neuron formula.")
        print("r = 1 / (tau_ref + tau_m * ln((mu - V_reset) / (mu - V_th)))")
        print("\nPlugging in the numbers from the problem:")
        print(f"r = 1 / ({tau_ref} + {tau_m} * ln(({mu} - {V_reset}) / ({mu} - {V_th})))")
        print(f"r = 1 / ({tau_ref} + {tau_m} * ln({log_arg_numerator} / {log_arg_denominator}))")
        print(f"r = 1 / ({tau_ref} + {tau_m} * ln({log_argument}))")
        print(f"r = {firing_rate_hz:.2f} Hz\n")

        print("Step 3: Final Answer")
        print(f"The firing rate of a typical neuron, rounded to the nearest integer, is {final_answer} Hz.")
        
    else:
        final_answer = 0
        print("The total input potential 'mu' is not above the threshold, so the neuron does not fire.")

    return final_answer

if __name__ == '__main__':
    answer = calculate_neuron_firing_rate()
    print(f"\n<<< {answer} >>>")