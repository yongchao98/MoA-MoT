import sys

# A full classification of homogeneous planar continua reveals a small, finite number of such spaces.
# This script lists these classes and counts them.

# 1. Define the known homeomorphism classes of homogeneous planar continua.
#    These are the Point, the Circle, the Pseudo-arc, and the Circle of pseudo-arcs.
#    Each represents one unique class.
classes = {
    "The Point": 1,
    "The Circle": 1,
    "The Pseudo-arc": 1,
    "The Circle of pseudo-arcs": 1
}

# 2. Print the name of each class.
print("The distinct homeomorphism classes of homogeneous planar continua are:")
for name in classes:
    print(f"- {name}")

# 3. To find the total number, we can sum the counts (1 for each class).
total_classes = sum(classes.values())

# 4. Display the calculation as a sum, showing each number in the final equation.
#    Create a string like "1 + 1 + 1 + 1".
sum_expression = " + ".join([str(v) for v in classes.values()])

print("\nThe total number is the sum of one for each class:")
print(f"{sum_expression} = {total_classes}")
