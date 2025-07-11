import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall device
    with terminals 3 and 4 floated using the Landauer-Büttiker formalism.
    """
    # Define symbols for the voltages and the fundamental constants
    # V is the source voltage at terminal 1.
    # V3 and V4 are the unknown voltages at the floating terminals.
    # G0 represents the quantum of conductance, e^2/h.
    V, V3, V4, G0 = sympy.symbols('V V3 V4 G0')

    # --- Landauer-Büttiker Setup ---
    # In a QSH system, there are two counter-propagating edge states (one for each spin)
    # between any two adjacent terminals. Let's label terminals 1, 2, 3, 4 clockwise.
    #
    # N_i: Number of channels leaving terminal i.
    # For each terminal, one channel leaves clockwise and one counter-clockwise.
    # N1=2 (to 2, 4), N2=2 (to 3, 1), N3=2 (to 4, 2), N4=2 (to 1, 3)
    N = {1: 2, 2: 2, 3: 2, 4: 2}

    # T_ij: Transmission from terminal j to terminal i.
    # This is 1 if a channel exists, 0 otherwise.
    # T_ij = 1 for adjacent terminals, 0 for opposite terminals.
    # T21=1 (from 1 to 2), T41=1 (from 1 to 4)
    # T12=1 (from 2 to 1), T32=1 (from 2 to 3)
    # T23=1 (from 3 to 2), T43=1 (from 3 to 4)
    # T14=1 (from 4 to 1), T34=1 (from 4 to 3)
    # T31, T13, T42, T24 are all 0.
    
    # --- Boundary Conditions ---
    # Terminal 1 is source: V1 = V
    # Terminal 2 is drain: V2 = 0
    # Terminals 3 and 4 are floated: I3 = 0, I4 = 0
    V1 = V
    V2 = 0

    # --- Landauer-Büttiker Equations for floated terminals ---
    # I_i = G0 * [ (sum_{j!=i} T_ij * V_j) - N_i * V_i ] = 0
    # For I3 = 0:
    # We have channels from terminal 2 to 3 (T32=1) and from 4 to 3 (T43=1)
    # N3 = 2
    # Eq3: G0 * [ (T31*V1 + T32*V2 + T34*V4) - N[3]*V3 ] = 0
    #      G0 * [ (0*V1 + 1*V2 + 1*V4) - 2*V3 ] = 0
    #      (V2 + V4) - 2*V3 = 0
    #      (0 + V4) - 2*V3 = 0
    eq3 = sympy.Eq(V4 - 2 * V3, 0)
    
    # For I4 = 0:
    # We have channels from terminal 1 to 4 (T41=1) and from 3 to 4 (T43=1)
    # N4 = 2
    # Eq4: G0 * [ (T41*V1 + T42*V2 + T43*V3) - N[4]*V4 ] = 0
    #      G0 * [ (1*V1 + 0*V2 + 1*V3) - 2*V4 ] = 0
    #      (V1 + V3) - 2*V4 = 0
    eq4 = sympy.Eq(V1 + V3 - 2 * V4, 0)

    # --- Solve for unknown voltages V3 and V4 ---
    # We have a system of two linear equations with two variables (V3, V4).
    # The third variable 'V' is treated as a known constant.
    solution = sympy.solve([eq3, eq4], (V3, V4))
    
    v3_val = solution[V3]
    v4_val = solution[V4]

    # --- Calculate the current I2 and the conductance G12 ---
    # The current at the drain (terminal 2) is:
    # I2 = G0 * [ (T21*V1 + T23*V3 + T24*V4) - N[2]*V2 ]
    #    = G0 * [ (1*V1 + 1*V3 + 0*V4) - 2*V2 ]
    #    = G0 * [ V1 + V3 - 2*V2 ]
    # Substitute V1=V, V2=0, and the solved value for V3:
    I2 = G0 * (V1 + v3_val - 2 * V2)

    # The two-terminal conductance G12 is defined as I2 / V1 (since V2=0).
    G12 = I2 / V1
    
    # Extract the numerical coefficient from the symbolic expression
    conductance_coeff = G12.subs(G0, 1)
    num, den = sympy.fraction(conductance_coeff)
    
    # Print the final result in a descriptive format
    print("The two-terminal conductance G12 is calculated as follows:")
    print(f"Voltage at floated terminal 3, V3 = {v3_val}")
    print(f"Voltage at floated terminal 4, V4 = {v4_val}")
    print(f"Current at terminal 2, I2 = G0 * (V + V3) = G0 * (V + {v3_val}) = {I2}")
    print("\nThe conductance G12 = I2 / V.")
    print(f"G12 = ({I2}) / V = {G12}")
    print("\nFinal Result:")
    # The prompt requests to output each number in the final equation.
    print(f"G_12 = ({num}/{den}) * e^2/h")

calculate_qsh_conductance()