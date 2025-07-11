def predict_energy_spectrum():
    """
    This function prints the energy spectrum for a harmonic oscillator with a
    quartic perturbation, as predicted by the one-loop self-energy diagram.
    The result is presented as a symbolic formula.
    """

    # --- Symbolic representation of the variables ---
    n = "n"
    h_bar = "ħ"
    omega_0 = "ω₀"
    u = "u"
    m = "m"
    En = f"E_{n}"

    # --- Breaking down the final equation into its parts ---
    # The structure is: E_n = part_A * part_B * sqrt(part_C)

    # Part A: The quantum number dependent term
    part_A = f"({n} + 1/2)"

    # Part B: Planck's constant
    part_B = h_bar

    # Part C: The renormalized frequency term (inside the square root)
    # It has the form: ω₀² + shift
    original_freq_sq = f"{omega_0}**2"
    shift_numerator = f"{u} * {h_bar}"
    shift_denominator = f"4 * {m}**2 * {omega_0}"
    freq_shift = f"({shift_numerator}) / ({shift_denominator})"
    
    # --- Assembling and printing the final equation ---
    print("The predicted energy spectrum, based on the self-energy diagram, is E_n = (n + 1/2)ħω', where ω' is the new renormalized frequency.")
    print("The final equation for the energy levels E_n is constructed as follows:\n")
    
    final_equation = f"{En} = {part_A} * {part_B} * sqrt({original_freq_sq} + {freq_shift})"
    
    print(f"Final Equation:")
    print(f"{final_equation}\n")

    print("Where each part of the equation represents:")
    print(f"1. The state-dependent part: {part_A}")
    print(f"2. Planck's constant: {part_B}")
    print(f"3. The squared original frequency: {original_freq_sq}")
    print(f"4. The numerator of the frequency shift: {shift_numerator}")
    print(f"5. The denominator of the frequency shift: {shift_denominator}")


predict_energy_spectrum()
