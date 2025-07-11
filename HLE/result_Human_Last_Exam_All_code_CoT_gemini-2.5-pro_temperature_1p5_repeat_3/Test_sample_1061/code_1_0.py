import sys

def model_tree_ring_isotope_change():
    """
    This function conceptually models the effect of the SE Asia monsoon
    on the 13C isotope ratio in Chinese pine trees.
    """
    # δ13C (delta-13C) is the standard notation for the 13C ratio, typically
    # expressed in per mil (‰). A more negative value means a lower 13C ratio.
    # Let's assume a representative initial δ13C value for the start of the period.
    initial_delta_13C = -25.0

    # The problem describes a declining 13C ratio, which points to a long-term
    # trend of increased water availability. This is driven by the SE Asia monsoon.
    # A stronger monsoon leads to a more negative δ13C value. We represent
    # this environmental effect as a negative change.
    # Note: This is a conceptual value for illustration.
    change_caused_by_monsoon_trend = -1.2

    # The final δ13C value is the sum of the initial state and the long-term change.
    final_delta_13C = initial_delta_13C + change_caused_by_monsoon_trend

    print("Explanation: A declining 13C ratio (more negative δ13C value) in tree rings indicates improved water availability, allowing the tree to discriminate more against the heavy 13C isotope.")
    print("The SE Asia monsoon is the main driver of water availability in this region.")
    print("\nConceptual Equation of the Predominant Factor:")
    print("Final_Ratio = Initial_Ratio + Change_From_Monsoon")

    # The following line prints out each number in the final conceptual equation,
    # demonstrating the decline.
    sys.stdout.write(f"{final_delta_13C:.1f}‰ = {initial_delta_13C:.1f}‰ + ({change_caused_by_monsoon_trend:.1f}‰)\n")

    print("\nAs the equation shows, a change driven by the monsoon results in a final ratio that has declined, matching the observations.")

model_tree_ring_isotope_change()