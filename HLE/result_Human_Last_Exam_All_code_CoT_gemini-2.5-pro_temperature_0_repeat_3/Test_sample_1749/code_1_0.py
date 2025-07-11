def analyze_abr_for_auditory_neuropathy():
    """
    This script simulates the logic used to identify a key feature of
    Auditory Neuropathy Spectrum Disorder (ANSD) on an ABR test.
    """

    # --- Simulated ABR Findings ---

    # In ANSD, a pre-neural response from the cochlea is present.
    # Let's represent the wave seen with a rarefaction stimulus.
    rarefaction_response = "Waveform A"

    # A key test is to see how the wave responds to an inverted stimulus (condensation).
    # For a Cochlear Microphonic (CM), the wave inverts.
    condensation_response_is_mirror_image = True

    # The duration of this inverting wave is also important.
    # A true CM is sustained, unlike a brief stimulus artifact.
    inverting_wave_duration_ms = 1.5 # A value greater than 1ms

    # In ANSD, the neural responses (ABR Waves I, III, V) are absent or abnormal.
    neural_waves_present = False

    print("Analyzing simulated ABR data for signs of Auditory Neuropathy...")
    print("-" * 50)

    # The diagnostic "equation" for the classic ANSD finding on ABR involves
    # checking for an inverting, long-duration wave in the absence of neural responses.
    if condensation_response_is_mirror_image and inverting_wave_duration_ms > 1.0 and not neural_waves_present:
        print("Result: Findings are consistent with Auditory Neuropathy.")
        print("\nExplanation of the key diagnostic marker (the Cochlear Microphonic):")
        # This section prints the components of the "final equation" or logic.
        print(f"1. The wave pattern in condensation is a mirror image of the wave pattern in rarefaction: {condensation_response_is_mirror_image}")
        print(f"2. The duration of this inverting pattern is > 1ms: {inverting_wave_duration_ms} ms")
        print(f"3. The neural ABR waves (I, III, V) are absent: {not neural_waves_present}")
        print("\nThis combination of findings points to a present Cochlear Microphonic with absent neural response, a hallmark of ANSD.")
    else:
        print("Result: Findings are not characteristic of the classic ANSD presentation.")

analyze_abr_for_auditory_neuropathy()