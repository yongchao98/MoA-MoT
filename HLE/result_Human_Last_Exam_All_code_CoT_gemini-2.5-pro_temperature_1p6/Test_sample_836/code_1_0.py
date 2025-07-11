def find_transform_name():
    """
    This function explains and identifies the space-time, double Fourier transform
    of the generalized pair correlation function in nuclear physics.
    """

    # Step 1: Define the function of interest conceptually.
    print("Step 1: The Generalized Pair Correlation Function")
    print("-------------------------------------------------")
    print("Let G(r, t) be the generalized pair correlation function.")
    print("This function measures the correlation between neutron populations at points")
    print("separated by a spatial vector 'r' and a time interval 't'.\n")

    # Step 2: Define the mathematical transformation.
    print("Step 2: The Space-Time, Double Fourier Transform")
    print("-------------------------------------------------")
    print("The goal is to find the name of the function S(k, omega) resulting")
    print("from the double Fourier transform of G(r, t).")
    print("The transformation is defined by the following integral:\n")
    
    # Printing the equation symbol by symbol as requested.
    equation_symbols = [
        "S(k, omega)", " = ", "Integral_over_t", "Integral_over_r",
        " G(r, t)", " * ", "exp[-i * (k.r - omega*t)]", " dr", " dt"
    ]
    print("Final Equation:")
    print(" ".join(equation_symbols))
    print("\nWhere:")
    print("  'k' is the wavevector (spatial frequency)")
    print("  'omega' is the angular frequency (temporal frequency)")
    print("  'i' is the imaginary unit\n")
    
    # Step 3: Provide the common name for the resulting function.
    print("Step 3: The Common Name")
    print("--------------------------")
    print("The resulting function, S(k, omega), is a spectral density that describes")
    print("the collective spatio-temporal dynamics of the neutron population.")
    print("In the nuclear criticality community, as in condensed matter physics,")
    print("this function is commonly called the:\n")
    
    final_answer = "Dynamic Structure Factor"
    
    print(f"*** {final_answer} ***")

if __name__ == "__main__":
    find_transform_name()