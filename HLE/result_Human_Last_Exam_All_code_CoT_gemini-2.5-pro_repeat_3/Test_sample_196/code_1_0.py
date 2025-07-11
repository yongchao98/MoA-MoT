import sympy
from sympy import pi, sqrt, integrate

# Define the symbolic variable for integration
r = sympy.symbols('r')

# Define the integration limit, which is the radius of the intersection circle
r_limit = sympy.Rational(3, 2)

# The volume integral is 2*pi times the integral over r.
# The integrand in r is: (y_cone - y_ellipsoid) * r
# V = 2 * pi * integral from 0 to 3/2 of [ (4 - 2r - 2*sqrt(1 - r^2/3)) * r ] dr
# V = 2 * pi * integral from 0 to 3/2 of [ 4r - 2r^2 - 2r*sqrt(1 - r^2/3) ] dr

# We calculate the integral by breaking it into three parts
# Part 1: Integral of 4r
integrand1 = 4 * r
integral_part1 = integrate(integrand1, (r, 0, r_limit))

# Part 2: Integral of -2r^2
integrand2 = -2 * r**2
integral_part2 = integrate(integrand2, (r, 0, r_limit))

# Part 3: Integral of -2r*sqrt(1 - r^2/3)
integrand3 = -2 * r * sqrt(1 - r**2 / 3)
integral_part3 = integrate(integrand3, (r, 0, r_limit))

# The total value of the integral over r is the sum of the parts
total_integral_r = integral_part1 + integral_part2 + integral_part3

# The final volume is 2*pi times this result
final_volume = 2 * pi * total_integral_r

# Print the breakdown of the final calculation
print("The volume calculation is based on the equation:")
print(f"Volume = 2 * pi * (Integral[Part 1] + Integral[Part 2] + Integral[Part 3])")
print("\nWhere the integrals are evaluated from r=0 to r=3/2:")
print(f"Integral[Part 1] = Integral[4r] dr = {integral_part1}")
print(f"Integral[Part 2] = Integral[-2r^2] dr = {integral_part2}")
print(f"Integral[Part 3] = Integral[-2r*sqrt(1 - r^2/3)] dr = {integral_part3}")

print("\nThe final equation with the calculated values is:")
print(f"Volume = 2 * pi * ( {integral_part1} + ({integral_part2}) + ({integral_part3}) )")
print(f"Volume = 2 * pi * ( {total_integral_r} )")
print(f"Volume = {final_volume}")