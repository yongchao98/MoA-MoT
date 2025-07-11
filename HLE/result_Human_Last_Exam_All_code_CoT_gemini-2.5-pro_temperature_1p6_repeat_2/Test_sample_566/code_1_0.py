def explain_broadband_cars():
    """
    Explains the principles of broadband CARS microscopy to determine the correct answer.
    """
    
    print("Thinking Process:")
    print("1. Coherent Anti-Stokes Raman Scattering (CARS) is a four-wave mixing process used for chemical imaging.")
    print("   - It uses pump (w_p) and Stokes (w_s) beams to excite a molecular vibration (Omega) when w_p - w_s = Omega.")
    print("   - A probe beam (w_pr, often w_p) then scatters off this vibration to generate a signal.")

    print("\n2. The signal generated is at the anti-Stokes frequency, w_as.")
    print("   - The key relationship is: w_as = w_p - w_s + w_pr.")
    print("   - Since Omega = w_p - w_s, and if w_pr = w_p, this simplifies to w_as = w_p + Omega.")

    print("\n3. In 'Broadband' CARS, one of the incident beams (e.g., the pump or Stokes) has a wide range of frequencies.")
    print("   - This allows the system to simultaneously match and excite an entire range of different vibrational frequencies (multiple Omega values) present in the sample.")

    print("\n4. Consequently, the generated anti-Stokes signal (w_as) is also broadband.")
    print("   - This broadband signal is a composite of different frequencies, where each frequency corresponds to a specific molecular vibration (Omega) that was excited.")

    print("\n5. This signal is then sent to a spectrometer, which separates the light into its constituent frequencies.")
    print("   - This allows for the measurement of a full vibrational spectrum in a single acquisition.")
    print("   - Therefore, the anti-Stokes beam contains 'distinguishable information' about the different vibrational modes.")

    print("\nConclusion:")
    print("The generated anti-Stokes beam is a spectrum that carries spectrally resolved information about the sample's molecular vibrations. This makes statement C the correct one.")

    print("\nOutputting the components of the final signal equation as requested:")
    # Printing each component of the equation w_as = w_p + Omega
    print("w_as")
    print("=")
    print("w_p")
    print("+")
    print("Omega")

explain_broadband_cars()

# The final answer is C.
print("<<<C>>>")