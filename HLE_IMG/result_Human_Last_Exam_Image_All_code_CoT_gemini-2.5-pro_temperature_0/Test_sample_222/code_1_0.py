import math
from fractions import Fraction

# The function is f(z) = z/(z-3/2) * Gamma(z), which can be rewritten as f(z) = Gamma(z+1)/(z-3/2).
# The singularities of f(z) are poles at z = 3/2 and at z = -1, -2, -3, ...
# The point z=0 is a removable singularity, so its residue is 0.
# The pole at z=3/2 is not enclosed by either contour.
# We need to calculate the residues at the poles z = -1, -2, -3.

# The residue of f(z) at z = -k (for k >= 1) is given by the formula:
# Res(f, -k) = (-2 / (2k + 3)) * ((-1)**(k-1) / (k-1)!)

# Calculate residues as exact fractions
res_minus_1 = Fraction(-2, 2*1 + 3) * Fraction((-1)**(1-1), math.factorial(1-1)) # k=1
res_minus_2 = Fraction(-2, 2*2 + 3) * Fraction((-1)**(2-1), math.factorial(2-1)) # k=2
res_minus_3 = Fraction(-2, 2*3 + 3) * Fraction((-1)**(3-1), math.factorial(3-1)) # k=3

# Winding numbers from the contours C1 and C2
# For C1: Ind_C1(-1)=+1, Ind_C1(-2)=-1, Ind_C1(-3)=+1
# For C2: Ind_C2(-1)=+1, Ind_C2(-2)=+1, Ind_C2(-3)=0
# Total winding numbers = Ind_C1 + Ind_C2
total_wind_minus_1 = 1 + 1
total_wind_minus_2 = -1 + 1
total_wind_minus_3 = 1 + 0

# The sum of the integrals is I = 2*pi*i * Sum[ TotalInd(zk) * Res(f, zk) ]
# The imaginary part is Im(I) = 2*pi * Sum[ TotalInd(zk) * Res(f, zk) ]

sum_of_terms = (total_wind_minus_1 * res_minus_1 +
                total_wind_minus_2 * res_minus_2 +
                total_wind_minus_3 * res_minus_3)

# Print the step-by-step calculation
print("The imaginary part of the sum of the integrals is calculated using the Residue Theorem.")
print("Im(I) = 2 * pi * Sum[ (Ind_C1(zk) + Ind_C2(zk)) * Res(f, zk) ]")
print("\nRelevant poles and their data:")
print(f"Pole z=-1: Residue = {res_minus_1}, Total Winding Number = {total_wind_minus_1}")
print(f"Pole z=-2: Residue = {res_minus_2}, Total Winding Number = {total_wind_minus_2}")
print(f"Pole z=-3: Residue = {res_minus_3}, Total Winding Number = {total_wind_minus_3}")

print("\nThe final equation for the imaginary part is:")
# Im(I) = 2 * pi * [ (2 * -2/5) + (0 * 2/7) + (1 * -1/9) ]
print(f"Im(I) = 2 * pi * [ ({total_wind_minus_1})*({res_minus_1}) + ({total_wind_minus_2})*({res_minus_2}) + ({total_wind_minus_3})*({res_minus_3}) ]")
term1 = total_wind_minus_1 * res_minus_1
term2 = total_wind_minus_2 * res_minus_2
term3 = total_wind_minus_3 * res_minus_3
print(f"Im(I) = 2 * pi * [ {term1} + {term2} + {term3} ]")
print(f"Im(I) = 2 * pi * ( {sum_of_terms} )")
final_fraction = 2 * sum_of_terms
print(f"Im(I) = ({final_fraction.numerator}/{final_fraction.denominator}) * pi")

# Calculate the final numerical value
numerical_value = float(final_fraction) * math.pi
print(f"\nThe numerical value is approximately: {numerical_value}")