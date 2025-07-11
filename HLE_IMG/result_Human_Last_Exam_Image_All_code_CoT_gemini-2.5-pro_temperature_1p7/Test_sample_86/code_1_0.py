import sympy

# Define symbols for mathematical representation
u, hbar, m, omega_0, n, x = sympy.symbols("u hbar m omega_0 n x")

# --- Introduction ---
print("This problem asks for the energy spectrum of a harmonic oscillator with a quartic perturbation,")
print("as predicted by the first-order self-energy diagram (tadpole diagram).")
print("Here is a step-by-step derivation leading to the final result.")
print("-" * 50)

# --- Step 1: Effective Potential Correction ---
print("\nStep 1: Determine the effective potential from the diagram.")
print("The interaction Hamiltonian is H_1 = (u / 4!) * x^4.")
print("The tadpole diagram corresponds to a mean-field approximation where two of the four x operators")
print("are replaced by their ground-state expectation value, <0|x^2|0>.")
print("The combinatorial factor for choosing 2 out of 4 operators is C(4,2) = 6.")
print("This gives an effective correction to the potential, delta_V:")
print("delta_V = (u / 24) * 6 * <0|x^2|0> * x^2 = (u / 4) * <0|x^2|0> * x^2")
print("-" * 50)

# --- Step 2: Ground State Expectation Value ---
print("\nStep 2: State the ground state expectation value of x^2.")
print("For an unperturbed harmonic oscillator, this value is a known result:")
print("<0|x^2|0> = hbar / (2 * m * omega_0)")
print("-" * 50)

# --- Step 3: Substitute and find delta_V ---
print("\nStep 3: Substitute the expectation value into delta_V.")
print("delta_V = (u / 4) * (hbar / (2 * m * omega_0)) * x^2")
print("delta_V = (u * hbar) / (8 * m * omega_0) * x^2")
print("-" * 50)

# --- Step 4: Find the new frequency ---
print("\nStep 4: Determine the new renormalized frequency.")
print("The new total potential V' is the sum of the original V_0 and the correction delta_V:")
print("V' = V_0 + delta_V = (1/2)*m*omega_0^2*x^2 + (u*hbar)/(8*m*omega_0)*x^2")
print("\nThis can be rewritten in the form of a new harmonic oscillator potential V' = (1/2)*m*omega^2*x^2:")
print("V' = [ (1/2)*m*omega_0^2 + (u*hbar)/(8*m*omega_0) ] * x^2")
print("By comparing the terms, we find the new frequency squared (omega^2):")
print("(1/2)*m*omega^2 = (1/2)*m*omega_0^2 + (u*hbar)/(8*m*omega_0)")
print("\nSolving for omega^2:")
print("omega^2 = omega_0^2 + (u * hbar) / (4 * m^2 * omega_0)")
print("-" * 50)

# --- Step 5: State the final energy spectrum ---
print("\nStep 5: Write the energy spectrum using the new frequency.")
print("The energy levels of a harmonic oscillator with frequency omega are E_n = (n + 1/2)*hbar*omega.")
print("Using our renormalized frequency, the predicted spectrum is:")
print("\nE_n = (n + 1/2) * hbar * sqrt( omega_0^2 + (u * hbar) / (4 * m^2 * omega_0) )")
print("-" * 50)

# --- Final Answer Output ---
print("\nThe final equation is:")
print("E_n", "=", "(", "n", "+", "1", "/", "2", ")", "*", "hbar", "*", "sqrt(", "omega_0", "^", "2", "+", "(", "u", "*", "hbar", ")", "/", "(", "4", "*", "m", "^", "2", "*", "omega_0", ")", ")")