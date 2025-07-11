import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a thin-walled, asymmetric channel section.

    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, it is located at a distance 'e'
    from the centerline of the web.

    This calculation assumes a thin-walled section where the thicknesses are small
    compared to the height and flange widths.
    """
    # Define the dimensions of the asymmetric channel section in mm
    h = 200.0   # Height of the web (between flange centerlines)
    tw = 5.0    # Thickness of the web
    b1 = 80.0   # Width of the top flange (from web centerline)
    t1 = 8.0    # Thickness of the top flange
    b2 = 100.0  # Width of the bottom flange (from web centerline)
    t2 = 10.0   # Thickness of the bottom flange

    print("Calculating the shear center location 'e' for an asymmetric channel section.")
    print("Dimensions (in mm):")
    print(f"  Web Height (h): {h}")
    print(f"  Web Thickness (tw): {tw}")
    print(f"  Top Flange Width (b1): {b1}")
    print(f"  Top Flange Thickness (t1): {t1}")
    print(f"  Bottom Flange Width (b2): {b2}")
    print(f"  Bottom Flange Thickness (t2): {t2}\n")

    # Formula for the shear center offset 'e' from the web centerline
    # e = (t1 * b1^2 + t2 * b2^2) / ( (2/3)*tw*h + 2*t1*b1 + 2*t2*b2 )

    # Calculate the numerator and denominator
    numerator = t1 * b1**2 + t2 * b2**2
    denominator = (2/3) * tw * h + 2 * t1 * b1 + 2 * t2 * b2

    # Calculate the shear center location 'e'
    e = numerator / denominator

    # Print the equation with the values substituted
    print("Formula: e = (t1*b1^2 + t2*b2^2) / ((2/3)*tw*h + 2*t1*b1 + 2*t2*b2)")
    print("\nSubstituting the values into the equation:")
    # Print the numerator calculation
    print(f"Numerator = ({t1} * {b1}**2) + ({t2} * {b2}**2) = {t1*b1**2} + {t2*b2**2} = {numerator}")
    # Print the denominator calculation
    print(f"Denominator = (2/3 * {tw} * {h}) + (2 * {t1} * {b1}) + (2 * {t2} * {b2}) = { (2/3) * tw * h :.2f} + {2 * t1 * b1} + {2 * t2 * b2} = {denominator:.2f}")

    print("\nFinal Equation:")
    print(f"e = {numerator} / {denominator:.2f}")

    # Print the final result
    print(f"\nThe shear center is located at a distance e = {e:.2f} mm from the centerline of the web.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()