import math

def calculate_cp_length():
    """
    Calculates the length of CP based on the given length of AB.
    The problem's geometric construction leads to the finding that
    the length of CP is the length of AB multiplied by the golden ratio.
    """
    # The given length of the segment AB.
    ab_length = 18468

    # The golden ratio, phi, is (1 + sqrt(5)) / 2.
    # The derived formula is: CP = AB * phi.
    cp_length = ab_length * ( (1 + math.sqrt(5)) / 2 )

    # As requested, here is the final equation with each number.
    # We display the formula with the input value of AB and the final calculated result for CP.
    print(f"{ab_length} * ( (1 + sqrt(5)) / 2 ) = {cp_length:.4f}")

    # The final answer in the special format <<<>>>
    print(f"<<<{cp_length:.4f}>>>")

calculate_cp_length()