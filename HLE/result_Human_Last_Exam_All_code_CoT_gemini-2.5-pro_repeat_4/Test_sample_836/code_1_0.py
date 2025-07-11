def main():
    """
    This function provides the name for the space-time, double Fourier transform
    of the generalized pair correlation function in nuclear criticality.
    It also presents the common symbolic equation for this quantity.
    """
    
    name = "Power Spectral Density (PSD)"
    
    print(f"The space-time, double Fourier transform of the generalized pair correlation function is commonly called the '{name}' of the neutron noise.")
    print("\n-------------------------------------------------")
    print("This quantity is central to reactor noise analysis.")
    print("In the simple case of a point reactor model (ignoring spatial effects), the Auto-Power Spectral Density (APSD) has a characteristic Lorentzian shape given by the following equation:")
    
    # The prompt asks to "output each number in the final equation".
    # As this is a symbolic formula, I will print and describe each component.
    
    print("\nEquation: APSD(w) = A / (alpha**2 + w**2)")
    
    print("\nWhere the components are:")
    print("  APSD(w): The Auto-Power Spectral Density as a function of angular frequency 'w'.")
    print("  A      : An amplitude or magnitude constant. It is not a single number but a factor dependent on the fission rate, detector efficiency, and the Diven factor (variance of neutrons per fission).")
    print("  alpha  : The prompt neutron decay constant, also known as the Rossi-alpha. This is a key parameter that characterizes the system's subcriticality.")
    print("  w      : The angular frequency (omega).")
    print("  **2    : Represents squaring the term.")

if __name__ == "__main__":
    main()
