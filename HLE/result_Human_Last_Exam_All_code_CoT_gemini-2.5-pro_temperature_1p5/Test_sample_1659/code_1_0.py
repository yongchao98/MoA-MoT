import numpy as np

def run_quantum_calculation():
    """
    Calculates the probability of measuring the state |100>
    from the final state |psi_4>.
    """
    # The basis states are ordered as |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
    # The final state from the derivation is:
    # |psi_4> = 0.5*|000> + 0.5*|100> + 0.5*|011> - 0.5*|111>
    psi_4 = np.array([0.5, 0, 0, 0.5, 0.5, 0, 0, -0.5])

    # The target outcome is |100>
    # Its index in the state vector is 1*2^2 + 0*2^1 + 0*2^0 = 4
    target_index = 4

    # Get the amplitude of the target state |100>
    amplitude = psi_4[target_index]

    # The probability is the squared magnitude of the amplitude
    probability = np.abs(amplitude)**2

    print("Final State |psi_4> = (1/2)|000> + (1/2)|100> + (1/2)|011> - (1/2)|111>")
    print(f"The amplitude of the target state |100> is: {amplitude}")
    print("The probability is the squared magnitude of the amplitude:")
    print(f"P(|100>) = |{amplitude}|^2 = {probability}")

run_quantum_calculation()