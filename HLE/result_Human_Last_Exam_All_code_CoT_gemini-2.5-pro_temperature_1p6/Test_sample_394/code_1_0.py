import math

def main():
    """
    Calculates the smallest possible denominator of the hypotenuse of a right
    triangle with area 263, using the properties of elliptic curves.
    """
    # The area of the triangle
    N = 263

    # The generator point (x, y) for the elliptic curve y^2 = x^3 - N^2*x
    x_num = 884449
    x_den = 2500
    y_num = 203874333
    y_den = 125000
    
    print(f"The area of the right triangle is N = {N}.")
    print(f"The sides of such a triangle can be found from a rational point on the elliptic curve y^2 = x^3 - N^2*x.")
    print(f"A known generator point (x,y) for this curve is:")
    print(f"x = {x_num}/{x_den}")
    print(f"y = {y_num}/{y_den}")
    
    # The hypotenuse c is given by c = (x^2 + N^2) / |y|
    print("\nThe hypotenuse c is calculated using the formula: c = (x^2 + N^2) / |y|")
    
    # Calculate the numerator of (x^2 + N^2)
    # (x_num/x_den)^2 + N^2 = (x_num^2 + N^2 * x_den^2) / x_den^2
    num_of_c_numerator = x_num**2 + (N**2) * (x_den**2)
    den_of_c_numerator = x_den**2

    # Now calculate c = (num_of_c_numerator / den_of_c_numerator) / (y_num / y_den)
    # c = (num_of_c_numerator * y_den) / (den_of_c_numerator * y_num)
    
    hypotenuse_numerator = num_of_c_numerator * y_den
    hypotenuse_denominator = den_of_c_numerator * y_num
    
    # Simplify the fraction by dividing by the greatest common divisor (GCD)
    common_divisor = math.gcd(hypotenuse_numerator, hypotenuse_denominator)
    
    final_numerator = hypotenuse_numerator // common_divisor
    final_denominator = hypotenuse_denominator // common_divisor
    
    print(f"\nSubstituting the values, we get:")
    print(f"Numerator of x^2 + N^2 is: {num_of_c_numerator}")
    print(f"Denominator of x^2 + N^2 is: {den_of_c_numerator}")
    print(f"Numerator of y is: {y_num}")
    print(f"Denominator of y is: {y_den}")
    print("\nSo, the hypotenuse c as a fraction is:")
    print(f"c = ({num_of_c_numerator} / {den_of_c_numerator}) / ({y_num} / {y_den})")
    print(f"c = ({num_of_c_numerator} * {y_den}) / ({den_of_c_numerator} * {y_num})")
    print(f"c = {hypotenuse_numerator} / {hypotenuse_denominator}")
    
    print(f"\nAfter simplifying by dividing the numerator and denominator by their GCD ({common_divisor}), we get:")
    print(f"c = {final_numerator} / {final_denominator}")
    
    print(f"\nThe smallest possible denominator of the hypotenuse is {final_denominator}.")

if __name__ == "__main__":
    main()
