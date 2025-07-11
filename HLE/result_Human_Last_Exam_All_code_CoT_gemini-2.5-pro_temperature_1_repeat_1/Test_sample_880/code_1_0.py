import pywt

def analyze_wavelet_candidates():
    """
    Analyzes candidate mother wavelets for rainfall time series analysis
    and recommends the best fit.
    """
    print("Analyzing properties of candidate mother wavelets...\n")

    # List of candidate wavelet short names from the answer choices
    candidates = {
        "A": "db1",   # Daubechies1 (Haar)
        "B": "sym2",  # Symlet2
        "C": "coif1", # Coiflet1
        "E": "db2"    # Daubechies2
    }

    # Signal characteristics
    signal_description = "Daily rainfall data is characterized by sharp, abrupt, and discontinuous spikes."
    
    # Store properties for comparison
    properties = {}

    for key, name in candidates.items():
        try:
            wavelet = pywt.Wavelet(name)
            properties[key] = {
                "name": wavelet.name,
                "family_name": wavelet.family_name,
                "support_length": wavelet.dec_len, # also known as filter length
                "symmetric": wavelet.symmetric,
                "orthogonal": wavelet.orthogonal
            }
        except ValueError:
            print(f"Could not retrieve info for wavelet: {name}")

    # Print the properties of each wavelet
    for key, prop in properties.items():
        print(f"--- Option {key} ({prop['name']}) ---")
        print(f"Family: {prop['family_name']}")
        print(f"Support Length: {prop['support_length']}")
        print(f"Symmetric: {prop['symmetric']}")
        print("-" * 25)

    print("\n--- Recommendation ---")
    print(f"Signal Type: {signal_description}")
    print("\nReasoning:")
    print("1. For signals with sharp transients like rainfall, excellent time localization is critical. This is achieved with a wavelet that has a small 'support length'.")
    print("2. The shape of the wavelet should ideally match the shape of the signal features. Rainfall spikes are similar to step functions.")
    
    print("\nConclusion:")
    print("The Daubechies1 ('db1') or Haar wavelet has the smallest possible support length (2), providing the best time localization to pinpoint rainfall events.")
    print("Its square, step-like shape is a great match for the abrupt nature of daily rainfall.")
    print("Smoother wavelets (sym2, coif1, db2) have longer support, which can blur the sharp peaks in time and are less suitable for this specific signal type.")

    final_choice = "A"
    print(f"\nTherefore, the best fit among the choices is Daubechies1.")
    
    # Final answer format as requested
    print(f"\n<<<{final_choice}>>>")

if __name__ == '__main__':
    analyze_wavelet_candidates()