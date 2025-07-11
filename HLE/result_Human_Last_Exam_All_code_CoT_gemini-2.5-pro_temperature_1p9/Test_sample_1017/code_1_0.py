from fractions import Fraction
import math

def main():
    """
    Computes the stable commutator length of the element g1*h2.
    """
    # Step 1: Define the rotation numbers from the problem description.
    g_rot = Fraction(2, 27)
    h_rot = Fraction(16, 27)

    # Step 2: State the formula for scl in the free product.
    print("The stable commutator length (scl) is computed using the formula:")
    print("scl(g_1 * h_2) = (1/2) * ||rot(g) - rot(h)||_Z")
    print("where ||x||_Z is the distance from x to the nearest integer.")
    print("")

    # Step 3: Calculate the difference in rotation numbers.
    delta_rot = g_rot - h_rot
    print(f"First, we calculate the difference of the rotation numbers:")
    print(f"rot(g) - rot(h) = {g_rot.numerator}/{g_rot.denominator} - {h_rot.numerator}/{h_rot.denominator} = {delta_rot.numerator}/{delta_rot.denominator}")
    print("")

    # Step 4: Find the distance to the nearest integer.
    print(f"Next, we find the distance from {delta_rot.numerator}/{delta_rot.denominator} to the nearest integer:")
    # A simple way to find the nearest integer m is to calculate floor(x + 0.5).
    # This works for both positive and negative x.
    nearest_integer = math.floor(float(delta_rot) + 0.5)
    
    # distance = abs(delta_rot - nearest_integer)
    dist_to_0 = abs(delta_rot)
    dist_to_neg_1 = abs(delta_rot - Fraction(-1))
    distance = min(dist_to_0, dist_to_neg_1)
    
    print(f"The nearest integers to {float(delta_rot):.4f} are 0 and -1.")
    print(f"Distance to 0: |{delta_rot.numerator}/{delta_rot.denominator}| = {dist_to_0.numerator}/{dist_to_0.denominator}")
    print(f"Distance to -1: |{delta_rot.numerator}/{delta_rot.denominator} - (-1)| = |{dist_to_neg_1.numerator}/{dist_to_neg_1.denominator}| = {dist_to_neg_1.numerator}/{dist_to_neg_1.denominator}")
    print(f"The minimum distance is {distance.numerator}/{distance.denominator}.")
    print("")

    # Step 5: Calculate the final scl value.
    scl = distance / 2
    print("Finally, the stable commutator length is half of this distance:")
    print(f"scl = (1/2) * {distance.numerator}/{distance.denominator} = {scl.numerator}/{scl.denominator}")

if __name__ == "__main__":
    main()
