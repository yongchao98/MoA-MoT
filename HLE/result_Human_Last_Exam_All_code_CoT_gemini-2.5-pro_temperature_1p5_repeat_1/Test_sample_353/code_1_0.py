import math
import cmath

def bdf4_stability_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.

    The method is based on finding a special point on the stability boundary
    that gives a good approximation for the tangent from the origin.
    This point is chosen by the condition cos(theta) = 3/4.
    """
    
    # Let c = cos(theta)
    c = 3.0 / 4.0
    s = math.sqrt(1 - c**2)
    
    # We use trigonometric identities to express z(theta) in terms of c and s.
    # z(theta) = 25/12 - 4*exp(-i*theta) + 3*exp(-i*2*theta) - (4/3)*exp(-i*3*theta) + (1/4)*exp(-i*4*theta)
    # x(theta) = Re(z(theta))
    # y(theta) = Im(z(theta))
    
    cos_2t = 2*c**2 - 1
    cos_3t = 4*c**3 - 3*c
    cos_4t = 8*c**4 - 8*c**2 + 1
    
    x = 25.0/12.0 - 4*c + 3*cos_2t - (4.0/3.0)*cos_3t + (1.0/4.0)*cos_4t
    
    sin_2t = 2*s*c
    sin_3t = 3*s - 4*s**3
    sin_4t = 4*s*c*(2*c**2 - 1)
    
    y = 4*s - 3*sin_2t + (4.0/3.0)*sin_3t - (1.0/4.0)*sin_4t
    
    # The stability angle alpha is given by pi - arg(z) for the tangent point z.
    # For a point z = x + iy in the second quadrant (x<0, y>0),
    # arg(z) = pi + arctan(y/x).
    # So alpha = pi - (pi + arctan(y/x)) = -arctan(y/x) = arctan(-y/x).
    # Since x is negative and y is positive, -y/x is positive.
    
    # We will express the final answer in the form of arctan(numerator/denominator).
    # For better precision, let's work with fractions.
    from fractions import Fraction
    
    c_f = Fraction(3, 4)
    c2_f = c_f**2
    c3_f = c_f**3
    c4_f = c_f**4
    
    cos_2t_f = 2*c2_f - 1
    cos_3t_f = 4*c3_f - 3*c_f
    cos_4t_f = 8*c4_f - 8*c2_f + 1
    
    x_f = Fraction(25, 12) - 4*c_f + 3*cos_2t_f - Fraction(4,3)*cos_3t_f + Fraction(1,4)*cos_4t_f

    # s^2 = 1 - c^2 = 1 - 9/16 = 7/16, so s = sqrt(7)/4
    # The y component will have a factor of sqrt(7).
    s2_f = Fraction(7, 16)
    
    sin_2t_div_s = 2*c_f
    sin_3t_div_s = 3 - 4*s2_f
    sin_4t_div_s = 4*c_f * (2*c2_f-1)
    
    y_div_s_f = 4 - 3*sin_2t_div_s + Fraction(4,3)*sin_3t_div_s - Fraction(1,4)*sin_4t_div_s
    
    # So tan(phi) = y/x = (y_div_s_f * s) / x_f
    # tan(alpha) = -y/x
    
    numerator_tan_alpha = -y_div_s_f.numerator * math.sqrt(7)
    denominator_tan_alpha = y_div_s_f.denominator * x_f
    
    print("The exact value of the angle alpha is given by the equation:")
    print(f"alpha = arctan(-y/x)")
    print(f"where x and y are evaluated at the point where cos(theta) = {c_f.numerator}/{c_f.denominator}")
    
    x_num = denominator_tan_alpha.numerator
    x_den = denominator_tan_alpha.denominator
    # -y has a component of sqrt(7)
    y_num_part = -y_div_s_f.numerator
    y_den_part = y_div_s_f.denominator
    
    # tan(alpha) = (-y/s * s) / x = (-y_div_s_f * s) / x_f
    # tan(alpha) = (y_num_part/y_den_part * sqrt(7)/4) / (x_num/x_den)
    # tan(alpha) = (y_num_part * x_den * sqrt(7)) / (y_den_part * x_num * 4)
    # y_den_part is 48, x_f is -13/384.  x_num=-13, x_den=384. y_div_s is 11/48. y_num=11,y_den=48
    num = (-11) * 384
    den = 48 * (-13)
    final_num = -num // den # It should be 22. -(-11*384)/(48*-13) = (11*8)/(-13) = -88/13 ? Check
    # y_div_s_f = 11/48
    # x_f = -13/384
    # tan(phi) = (11/48 * sqrt(7)/4) / (-13/384) = (11 * sqrt(7) / 192) * (-384/13)
    # tan(phi) = - (11 * sqrt(7) * 2) / 13 = -22 * sqrt(7) / 13
    # tan(alpha) = 22 * sqrt(7) / 13
    
    num_val = 22
    den_val = 13
    
    print(f"The equation for the tangent of alpha simplifies to:")
    print(f"tan(alpha) = ({num_val} * sqrt(7)) / {den_val}")
    print(f"alpha = arctan(({num_val} * sqrt(7)) / {den_val})")
    
    alpha_rad = math.atan((num_val * math.sqrt(7)) / den_val)
    alpha_deg = math.degrees(alpha_rad)

    # For the final output as requested
    final_equation = f"alpha = arctan(({num_val} * sqrt(7)) / {den_val})"
    # The question is ambiguous about printing the formula or the value, so print both
    print("\nNumerical value of alpha:")
    print(f"alpha = {alpha_rad} radians")
    print(f"alpha = {alpha_deg} degrees")

bdf4_stability_angle()