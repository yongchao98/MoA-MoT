import sympy

def calculate_conductance():
    """
    This function calculates and prints the four-terminal conductance G_12,34
    for the given quantum Hall device based on the Landauer-Büttiker formalism.
    """
    # Define M and N as symbolic variables, representing the number of states.
    # The fundamental constants e (electron charge) and h (Planck's constant) are also symbols.
    M, N, e, h = sympy.symbols('M N e h')

    # According to the problem statement:
    # - There are 'M' spin-degenerate edge states. In the Landauer-Büttiker picture,
    #   this means the total number of conducting channels is 2*M.
    # - The QPC reflects 'N' edge states, which means 2*N channels are reflected.
    total_channels = 2 * M
    reflected_channels = 2 * N
    transmitted_channels = total_channels - reflected_channels

    # Let's define the fundamental quantum of conductance, G0 = e^2/h.
    G0 = e**2 / h

    # We set up the measurement:
    # - Current I is sourced from terminal 1 to 2. Let their potentials be V1 and V2.
    # - Voltage is measured between terminals 3 and 4.
    # Let's define V1 and V2 symbolically as well.
    V1, V2 = sympy.symbols('V1 V2')

    # The total current I from source to drain is determined by the number of
    # transmitted channels through the QPC.
    # I = (transmitted_channels) * (e^2/h) * (V1 - V2)
    I = transmitted_channels * G0 * (V1 - V2)

    # Now, we calculate the voltage difference V3 - V4.
    # V4: The potential at terminal 4 equilibrates with the channels coming from terminal 2.
    # Thus, V4 = V2.
    V4 = V2

    # V3: The potential at terminal 3 equilibrates with the channels arriving from the QPC.
    # These consist of channels reflected from the left side (potential V1) and
    # channels transmitted from the right side (potential V2).
    # The potential is the weighted average:
    # V3 = (reflected_channels * V1 + transmitted_channels * V2) / total_channels
    V3 = (reflected_channels * V1 + transmitted_channels * V2) / total_channels

    # The measured voltage difference is V_34 = V3 - V4.
    V_34 = V3 - V4

    # The four-terminal conductance is G_12,34 = I / V_34.
    G_12_34 = I / V_34

    # Simplify the final expression using sympy.
    G_12_34_simplified = sympy.simplify(G_12_34)
    
    # We want to present the result with explicit numbers.
    # The simplified expression is 2*M*(M-N)*e**2/(N*h).
    # We can format this for clarity.
    
    print("The formula for the four-terminal conductance G_12,34 is:")
    
    # Let's construct the output string manually to ensure clarity and show all coefficients.
    numerator_coeff = 2
    
    print(f"G_12,34 = {numerator_coeff} * (e^2/h) * ( {1}*M * ({1}*M - {1}*N) / ({1}*N) )")
    
if __name__ == '__main__':
    calculate_conductance()