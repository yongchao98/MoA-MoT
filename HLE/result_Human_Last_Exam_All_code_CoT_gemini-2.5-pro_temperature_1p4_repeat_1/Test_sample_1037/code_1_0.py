import math

class Beth:
    """
    A class to represent Beth numbers for standardized string formatting.
    Beth(0) = aleph_0 (cardinality of natural numbers)
    Beth(1) = 2^aleph_0 (cardinality of the continuum)
    Beth(n) = 2^Beth(n-1)
    """
    def __init__(self, index):
        if not isinstance(index, int) or index < 0:
            raise ValueError("Beth index must be a non-negative integer.")
        self.index = index

    def __str__(self):
        return f"Beth_{self.index}"

# Step 1: Determine the cardinality of S/A
# A is the initial object in the category of scales, which is (Z, id). Its group is Z.
# S is the scale of the inclusion of Z into the hyperreals, *R. Its group is *R.
# The canonical map A -> S is the inclusion Z -> *R.
# The quotient S/A corresponds to the group quotient *R / Z.
# The cardinality of the hyperreals |*R| is |R^N| = (2^aleph_0)^aleph_0 = 2^aleph_0 = Beth_1.
# The cardinality of the integers |Z| is aleph_0 = Beth_0.
# For infinite groups G and H, |G/H| * |H| = |G|.
# So, |*R / Z| * |Z| = |*R|, which means |*R / Z| * Beth_0 = Beth_1.
# This implies |*R / Z| = Beth_1.
cardinality_S_div_A = Beth(1)

# Step 2: Determine the cardinality of B/S
# B is the terminal object. This corresponds to the trivial group {0}.
# S is the group *R.
# The canonical map S -> B is the unique zero map *R -> {0}. Its image is {0}.
# The quotient B/S corresponds to the group quotient {0} / Im(S->B) = {0} / {0}.
# The quotient group {0}/{0} has only one element.
cardinality_B_div_S = 1

# Step 3: Determine the cardinality of H_1(B/A, Q)
# The space B/A is the quotient of the underlying groups {0} / Z.
# The canonical map A -> B is the zero map Z -> {0}, so its image is {0}.
# The space B/A is {0} / {0}, which is a single-point space.
# We need the cardinality of the first homology group of a point, H_1(point, Q).
# For any n > 0, the homology group H_n(point, G) is the trivial group {0}.
# Therefore, H_1(B/A, Q) is {0}.
# The cardinality of the trivial group {0} is 1.
cardinality_H1 = 1

# Print the final result in the format "value1 value2 value3"
print(f"{cardinality_S_div_A} {cardinality_B_div_S} {cardinality_H1}")