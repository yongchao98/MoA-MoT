def calculate_qpc_conductance():
    """
    This function calculates and displays the formula for the four-terminal
    conductance G_12,34 of the specified quantum Hall device.

    Variables:
    M: The total number of spin-degenerate edge states.
    N: The number of edge states reflected by the quantum point contact (QPC).
    """

    # Symbolic representation of the variables
    M_symbol = "M"
    N_symbol = "N"

    print("Based on the Landauer-Büttiker formalism, the four-terminal conductance G_12,34 is derived.")
    print("The final formula is given in terms of M and N, and the quantum of conductance, e^2/h.")
    print("-" * 20)

    # The final derived formula for the conductance G_12,34
    # G = (I_12) / (V_34) = [ (M-N)*(e^2/h)*(V1-V2) ] / [ (N/M)*(V1-V2) ]
    # G = (M * (M - N) / N) * (e^2/h)

    # To satisfy the prompt "output each number in the final equation!",
    # we print the equation constructed with its symbolic components.
    print("Final Equation:")
    print(f"G_12,34 = ({M_symbol} * ({M_symbol} - {N_symbol})) / {N_symbol} * (e^2/h)")
    print("-" * 20)
    print("Where:")
    print(f"  '{M_symbol}' is the total number of edge states.")
    print(f"  '{N_symbol}' is the number of reflected edge states.")
    print("  'e^2/h' is the fundamental conductance quantum.")


if __name__ == '__main__':
    calculate_qpc_conductance()