def explain_fourier_transform_in_nuclear_criticality():
    """
    This function explains the space-time, double Fourier transform of the
    generalized pair correlation function as used in the nuclear criticality community.
    """

    # The common name for the transformed function
    common_name = "Power Spectral Density (PSD)"
    
    # A more specific name when considering two different locations/detectors
    specific_name = "Cross-Power Spectral Density (CPSD)"

    # The mathematical representation of the transformation
    # We describe the transform of the covariance function, which is directly
    # related to the pair correlation function.
    equation = "S(k, ω) = ∫∫ Cov(r, t) * e^(-i(k·r - ωt)) dr dt"

    print("The space-time, double Fourier transform of the generalized pair correlation function (or more precisely, the related covariance function) is commonly called the:\n")
    print(f"**{common_name}**\n")
    print(f"When specifically referring to the correlation between two different detectors, it is called the **{specific_name}**.\n")

    print("--- Mathematical Formulation ---")
    print("The final equation representing this relationship is:")
    print(f"  {equation}\n")

    print("--- Explanation of Symbols in the Final Equation ---")
    print("S(k, ω): The Power Spectral Density. This is the final transformed function.")
    print("      k: The wave vector, which is the variable in the spatial Fourier domain (units of 1/length).")
    print("      ω: The angular frequency, which is the variable in the temporal Fourier domain (units of 1/time).")
    print("\nCov(r, t): The covariance of the neutron flux fluctuations. It measures how the fluctuation at one point in space-time is related to another.")
    print("      r: The spatial separation vector between two points.")
    print("      t: The time lag between two points.")
    print("\ni: The imaginary unit, sqrt(-1).")
    print("\n∫∫ dr dt: Represents the double integral over all space (r) and all time (t).")
    print("\ne^(-i(k·r - ωt)): The complex exponential kernel of the Fourier transform.")

# Execute the function to print the explanation
explain_fourier_transform_in_nuclear_criticality()