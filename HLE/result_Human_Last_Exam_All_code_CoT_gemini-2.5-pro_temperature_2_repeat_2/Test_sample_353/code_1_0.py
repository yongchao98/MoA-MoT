import math
from fractions import Fraction

# Step 1: Calculate the coordinates x and y of the tangency point.
# Using the derived trigonometric values.
# cos_theta = 1/3
# cos_2theta = -7/9
# cos_3theta = -23/27
# cos_4theta = 17/81

x_t = (Fraction(25, 12) - 
       Fraction(4) * Fraction(1, 3) + 
       Fraction(3) * Fraction(-7, 9) - 
       Fraction(4, 3) * Fraction(-23, 27) + 
       Fraction(1, 4) * Fraction(17, 81))

# sin_theta = 2*sqrt(2)/3
# sin_2theta = 4*sqrt(2)/9
# sin_3theta = -10*sqrt(2)/27
# sin_4theta = -56*sqrt(2)/81

# Factor out sqrt(2) for y component
y_t_coeff = (Fraction(4) * Fraction(2, 3) - 
             Fraction(3) * Fraction(4, 9) + 
             Fraction(4, 3) * Fraction(-10, 27) - 
             Fraction(1, 4) * Fraction(-56, 81))

# x_t simplifies to -32/81
# y_t_coeff simplifies to 82/81
# So, y_t = (82*sqrt(2))/81

# Step 2: Calculate tan(alpha) = -y/x
# The sqrt(2) is handled symbolically.
# Note: we need to print each number in the final equation.
tan_alpha_num = -y_t_coeff.numerator * x_t.denominator
tan_alpha_den = y_t_coeff.denominator * x_t.numerator
# tan_alpha = -(82/81) / (-32/81) = 82/32 = 41/16
tan_alpha_frac = Fraction(tan_alpha_num, tan_alpha_den)

# The final angle is arctan( (41 * sqrt(2)) / 16 )
num1 = tan_alpha_frac.numerator
num2 = 2
num3 = tan_alpha_frac.denominator

print(f"The exact value for the angle alpha is arctan(({num1} * sqrt({num2})) / {num3})")
print("This corresponds to the equation:")
print(f"tan(alpha) = ({num1} * sqrt({num2})) / {num3}")
