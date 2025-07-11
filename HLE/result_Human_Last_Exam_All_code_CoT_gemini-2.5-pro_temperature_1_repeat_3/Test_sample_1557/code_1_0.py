def solve_wellington_mcq():
    """
    This function analyzes the historical statements about the Duke of Wellington's career
    and prints the numbers of the correct statements in the required format.
    """

    # Detailed analysis led to the identification of the following correct statements:
    #
    # Statement 1: Correct. Wellington's mastery of logistics in the Peninsular War was
    # built upon the principles of the mobile commissariat system he developed in India.
    # His ability to sustain an army in the field was a direct evolution of his Indian experience.
    #
    # Statement 6: Correct. The practice of integrating local auxiliary forces, which Wellington
    # used effectively in India, became a cornerstone of British imperial expansion. It was a
    # cost-effective and militarily sound strategy replicated across the empire in the 1830s-40s.
    #
    # Statement 8: Correct. The "flying column," a fast-moving, self-sufficient unit,
    # was a tactical innovation Wellington developed in India. This concept was successfully
    # adapted for both the large-scale Peninsular War and later colonial conflicts, such
    # as the First Anglo-Burmese War (1824-1826), where mobility was crucial.
    #
    # The other statements are incorrect for the following reasons:
    # - 2: His intelligence system was highly successful in Europe (Peninsular War).
    # - 3: Army-wide standardization did not occur in 1815; major reforms came much later.
    # - 4: Sandhurst was founded in 1801, not 1829, and its curriculum was Euro-centric.
    # - 5: The logistical systems were highly influential in other colonial campaigns.
    # - 7: His Peninsular commissariat was heavily influenced by his Indian experience.
    # - 9: The 1813 Charter Act was driven by economic and political, not military, principles.
    # - 10: His approach to civilian administration in India was a direct model for his actions in the Peninsula.

    correct_options = [1, 6, 8]

    # The numbers are already sorted, but we ensure it.
    correct_options.sort()

    # The final answer is the sorted numbers separated by a comma.
    final_answer = ",".join(map(str, correct_options))

    print(final_answer)

solve_wellington_mcq()