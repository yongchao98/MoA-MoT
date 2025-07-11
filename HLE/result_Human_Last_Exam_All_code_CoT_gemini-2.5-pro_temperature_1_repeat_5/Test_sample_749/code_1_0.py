import math

def calculate_hitting_probability():
    """
    Calculates the probability of a particle reaching site 0 in the limiting case (h->0).
    """
    # As h -> 0, all sites become blue with probability 1 for any finite region.
    # In a blue environment, the particle jumps left with probability 1/5
    # and right with probability 4/5.
    p_left = 1/5
    p_right = 4/5

    # The particle starts at site 3.
    start_site = 3

    # For a 1D random walk with a drift to the right (p_right > p_left),
    # the probability of reaching site 0 from a starting site k > 0 is (p_left / p_right)^k.
    if p_right <= p_left:
        # This case doesn't apply here but is included for completeness.
        # If drift is not to the right, hitting probability to the left is 1.
        probability = 1.0
    else:
        ratio = p_left / p_right
        probability = math.pow(ratio, start_site)

    print("This code calculates the probability that the initial particle ever reaches site 0")
    print("in the limiting environment (as h->0, all sites are blue).")
    print("This is an illustrative calculation, as the problem asks for the probability")
    print("of *infinite* visits, which is 0.")
    print("\n--- Equation ---")
    print(f"Let k be the starting site, k = {start_site}")
    print(f"Let p_L be the probability of jumping left, p_L = {p_left}")
    print(f"Let p_R be the probability of jumping right, p_R = {p_right}")
    print(f"P(hit 0 from k) = (p_L / p_R)^k")
    print(f"P(hit 0 from {start_site}) = ({p_left} / {p_right})^{start_site}")
    print(f"                   = ({p_left / p_right})^{start_site}")
    print(f"                   = {probability}")
    print("----------------")

calculate_hitting_probability()