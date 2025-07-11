import sympy

def calculate_conductance():
    """
    This function calculates and prints the four-terminal conductance G_12,34
    for the given quantum Hall device.
    """

    # Define the symbolic variables used in the problem.
    # M is the total number of spin-degenerate edge states.
    # N is the number of reflected edge states.
    # e is the elementary charge.
    # h is Planck's constant.
    M, N, e, h = sympy.symbols('M N e h')

    # The problem asks for the four-terminal conductance G_12,34.
    # Using the Landauer-Büttiker formalism, we derive the following relations:
    # 1. Voltage at probe 3: V_3 = (N/M)*V_1 + ((M-N)/M)*V_2
    # 2. Voltage at probe 4: V_4 = V_2
    # 3. Injected current: I = (M * e**2 / h) * ((M-N)/M) * (V_1 - V_2)

    # From these, we calculate the measured voltage difference V_3 - V_4:
    # V_3 - V_4 = ((N/M)*V_1 + ((M-N)/M)*V_2) - V_2
    #            = (N/M)*V_1 + ((M-N)/M - 1)*V_2
    #            = (N/M)*V_1 - (N/M)*V_2
    #            = (N/M) * (V_1 - V_2)

    # The conductance is G_12,34 = I / (V_3 - V_4).
    # G_12,34 = [ (M * e**2 / h) * ((M-N)/M) * (V_1 - V_2) ] / [ (N/M) * (V_1 - V_2) ]
    #         = [ (e**2 / h) * (M - N) ] / [ N/M ]
    #         = (M * (M - N) / N) * (e**2 / h)

    # We now print the final result.
    # The final equation gives the conductance G_12,34 in terms of M, N, and the conductance quantum e^2/h.
    # We will print out each component of the final formula as requested.
    
    m_str = "M"
    n_str = "N"
    
    # Numerator of the coefficient part of the equation
    numerator_coeff = f"1 * {m_str} * (1 * {m_str} - 1 * {n_str})"
    
    # Denominator of the coefficient part
    denominator_coeff = f"1 * {n_str}"
    
    # Fundamental constants part
    g0 = "e^2/h"

    print("The final equation for the four-terminal conductance G_12,34 is derived as:")
    print(f"G_12,34 = ( {numerator_coeff} ) / ( {denominator_coeff} ) * ( {g0} )")
    
    # Simpler visual form
    print("\nSimplified form:")
    print(f"G_12,34 = ({m_str}*({m_str} - {n_str})/{n_str}) * (e^2/h)")


if __name__ == "__main__":
    calculate_conductance()