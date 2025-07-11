import math

def find_computational_factor():
    """
    This function identifies the computational factor from historical CFD literature
    for the enthalpy-porosity method.
    """

    # The Carman-Kozeny source term in the enthalpy-porosity method is often modeled as:
    # S = -A_mush * ((1-beta)^2 / (beta^3 + c)) * u
    # where A_mush is the computational factor in question.

    # According to the foundational paper by Voller and Prakash (1987),
    # which was a simulation-only work, the value for the mushy region constant
    # (the computational factor) was specified.
    # This is the "prior published implementation" the question asks about.
    original_factor = 1.6 * math.pow(10, 6)

    # In the subsequent paper by Brent, Voller, and Reid (1988), which benchmarked
    # the method against the melting of gallium, this factor was modified to
    # better match experimental results.
    modified_factor = 1.6 * math.pow(10, 3)

    print("Identifying the original computational factor for the Carman-Kozeny source term in the enthalpy-porosity method.")
    print("-" * 60)
    print("The question asks for the value from the prior, simulation-only work, before it was adjusted for experimental benchmarks.")
    print("\nBased on a review of the relevant scientific literature:")
    print(f"1. The value used in the prior simulation-only paper (Voller & Prakash, 1987) was: {original_factor:.1e}")
    print(f"2. This value was later modified in the gallium melting benchmark paper (Brent, Voller & Reid, 1988) to: {modified_factor:.1e}")
    print("-" * 60)
    print("\nThe value originally used for this computational factor in the prior published implementation was therefore:")
    # We need to output each number in the final equation.
    # The value is 1.6 * 10^6.
    base = 1.6
    exponent = 6
    final_value = base * (10**exponent)
    
    print(f"{base} * 10^{exponent} = {final_value:.1e}")

find_computational_factor()