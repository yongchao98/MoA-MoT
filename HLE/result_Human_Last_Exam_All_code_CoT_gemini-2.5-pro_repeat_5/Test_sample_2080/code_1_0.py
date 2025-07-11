import sys

def calculate_steric_b22():
    """
    Calculates the second osmotic virial coefficient (B22) from steric-only
    behavior for a typical monoclonal antibody (mAb).
    """

    # The partial specific volume (v_bar) is a measure of the molecule's volume
    # in solution. For a typical mAb, this value is approximately 0.73 mL/g.
    partial_specific_volume = 0.73

    # In the hard-sphere approximation, the steric contribution to B22
    # is calculated as 4 times the partial specific volume. This represents
    # the repulsive force from the molecule's excluded volume.
    steric_factor = 4

    # The experimental B22 value of -7.585 mL/g results from the sum of this
    # positive steric term and a larger, negative term from attractive forces
    # (e.g., electrostatic and van der Waals interactions). We only need to
    # calculate the steric part.
    b22_steric = steric_factor * partial_specific_volume

    # Print the final equation showing all the numbers used, as requested.
    print("The calculation for the steric-only second osmotic virial coefficient (B22, steric) is:")
    print(f"B22, steric = {steric_factor} * {partial_specific_volume} mL/g")
    print(f"Result: {b22_steric} mL/g")

    # The final answer is returned for the wrapper.
    # The 'file=sys.stderr' argument is used to prevent this from being part of the standard output.
    print(f"\n<<<{b22_steric}>>>", file=sys.stderr)

if __name__ == '__main__':
    calculate_steric_b22()