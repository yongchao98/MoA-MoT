import math

def calculate_steric_b22():
    """
    Calculates the second osmotic virial coefficient from Steric-only behavior.

    The steric-only contribution (B22_steric) is calculated as four times
    the partial specific volume of the protein. A typical value for the
    partial specific volume of a monoclonal antibody is 0.73 mL/g.
    """
    # The partial specific volume for a typical protein/mAb in mL/g
    partial_specific_volume = 0.73

    # The factor relating partial specific volume to the steric B22
    factor = 4

    # Calculate the steric-only B22
    b22_steric = factor * partial_specific_volume

    print("The second osmotic virial coefficient from Steric-only behavior (B22_steric) is calculated from the partial specific volume (v_bar).")
    print("The relationship is: B22_steric = 4 * v_bar")
    print("\nUsing the standard value for a monoclonal antibody:")
    # Using round to ensure the output matches the expected precision from the calculation.
    print(f"B22_steric = {factor} * {partial_specific_volume} mL/g = {round(b22_steric, 3)} mL/g")

if __name__ == "__main__":
    calculate_steric_b22()