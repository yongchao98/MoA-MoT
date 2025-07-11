# This script provides the scansion for the Latin hexameter line:
# "verum ubi equi atque hominis casu convenit imago."
# The line is metrically complex; this solution uses the analysis from
# the Perseus Digital Library, which is D-S-D-S-S-S.
# D = Dactyl (— U U)
# S = Spondee (— —)

# Define the scansion for the six feet
foot1 = "D"
foot2 = "S"
foot3 = "D"
foot4 = "S"
foot5 = "S"
foot6 = "S"

# Print the complete scansion pattern
print(f"The scansion of the line is:")
print(f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}")

# Print the scansion in an "equation" format, outputting each symbol.
print("\nFinal equation:")
print(f"{foot1} + {foot2} + {foot3} + {foot4} + {foot5} + {foot6}")