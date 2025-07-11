from fractions import Fraction

# For an ensemble of disordered Majorana wires at the critical point, the
# transport properties are described by the BDI symmetry class. The dimensionless
# conductance 'g' follows the arcsine distribution.

# The statistical moments <g^n> can be calculated from this distribution.
# The average conductance <g> (the first moment, n=1) and the fourth moment <g^4> (n=4)
# have been derived to be the following exact fractions.

# Average value of the dimensionless conductance
avg_g = Fraction(1, 2)

# Fourth moment of the dimensionless conductance
fourth_moment_g = Fraction(35, 128)

# The problem asks for the ratio of the fourth moment to the average value.
ratio = fourth_moment_g / avg_g

print("Calculation of the ratio between the fourth and first statistical moments of conductance.")
print("-" * 80)
print(f"The average value (first moment) of the dimensionless conductance is: <g> = {avg_g}")
print(f"The fourth statistical moment of the dimensionless conductance is: <g^4> = {fourth_moment_g}")
print("-" * 80)
print("The final equation for the ratio is:")
print(f"Ratio = <g^4> / <g> = ({fourth_moment_g}) / ({avg_g}) = {ratio}")
print(f"\nThe numerical value of the ratio is: {float(ratio)}")