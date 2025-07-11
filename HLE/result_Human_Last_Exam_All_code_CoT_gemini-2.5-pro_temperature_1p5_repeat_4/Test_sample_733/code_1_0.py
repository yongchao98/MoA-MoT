# This problem is a known classification problem in group theory.
# The finite groups containing a maximal product-free set of size 2 have been fully classified.
# My task is to list these groups and provide a count.
# The groups are:
# 1. C4 (Cyclic group of order 4)
# 2. C5 (Cyclic group of order 5)
# 3. C6 (Cyclic group of order 6)
# 4. C2 x C2 (Klein four-group)
# 5. S3 (Symmetric group on 3 letters, also known as D3)
# 6. Q8 (Quaternion group)
#
# These 6 groups represent distinct isomorphism classes.

# The list of groups satisfying the condition
groups = [
    "C4 (Cyclic group of order 4)",
    "C5 (Cyclic group of order 5)",
    "C6 (Cyclic group of order 6)",
    "C2 x C2 (Klein four-group)",
    "S3 (Symmetric group of order 6)",
    "Q8 (Quaternion group)"
]

print("The finite groups which contain a maximal by inclusion product-free set of size 2 are:")
# This satisfies the requirement to "output each number in the final equation"
# by listing each component (group) being counted.
for i, group_name in enumerate(groups):
    print(f"{i+1}. {group_name}")

# The total number of such non-isomorphic groups.
count = len(groups)
print(f"\nThe total number of such groups is: {count}")