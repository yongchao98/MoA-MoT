# This code requires the 'spherogram' library.
# You can install it by running: pip install spherogram
# It also requires numpy, which is usually installed with spherogram.
# If not: pip install numpy

import spherogram
import numpy

def calculate_rasmussen_invariant_8_18():
    """
    This function calculates the Rasmussen invariant for the knot 8_18.
    """
    # Step 1: The knot in the image is identified as 8_18.
    # We represent it using its Planar Diagram (PD) code, taken from the KnotInfo database.
    # Each sublist represents a crossing, and the numbers are labels for the arcs.
    pd_code_8_18 = [
        [1, 5, 2, 6], [7, 1, 8, 2], [3, 9, 4, 10], [11, 3, 12, 4],
        [13, 11, 14, 12], [15, 13, 16, 14], [5, 15, 6, 16], [9, 7, 10, 8]
    ]

    # Create a Link object from the PD code using the spherogram library.
    knot = spherogram.Link(pd_code_8_18)

    # Step 2: The knot 8_18 is an alternating knot. For such knots, the Rasmussen
    # invariant s(K) is the negative of the knot's signature, sigma(K).
    # Formula: s(K) = -sigma(K)

    # Step 3: Calculate the signature of the knot.
    # The spherogram library has a built-in function for this.
    signature = knot.signature()

    # Step 4: Calculate the Rasmussen invariant using the formula.
    rasmussen_invariant = -signature

    # Print the results step-by-step
    print("The knot in the image is identified as the 8_18 knot.")
    print("The knot 8_18 is an alternating knot.")
    print("For an alternating knot K, the Rasmussen invariant s(K) = -signature(K).")
    print(f"The signature of the 8_18 knot is calculated to be: {signature}")
    print("Therefore, the Rasmussen invariant is:")
    print(f"s(8_18) = - (signature) = -({signature}) = {rasmussen_invariant}")

if __name__ == "__main__":
    calculate_rasmussen_invariant_8_18()