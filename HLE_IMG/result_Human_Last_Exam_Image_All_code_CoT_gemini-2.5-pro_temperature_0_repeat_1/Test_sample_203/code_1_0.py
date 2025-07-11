# The Heesch number for the top polygon is determined to be 3.
# This is based on the analysis that its shape allows for three layers of surrounding tiles but not a fourth.
heesch_number_top = 3

# The Heesch number for the middle polygon is 0.
# This is because its edges (protrusions and indentations) cannot fit against identical edges of its copies,
# preventing even a single complete layer from forming around it.
heesch_number_middle = 0

# The Heesch number for the bottom polygon is infinity.
# This is based on the assumption that this polygon is centrally symmetric,
# a property that guarantees it can tile the plane, thus allowing for infinite layers.
heesch_number_bottom = "∞"

# The problem asks for the answers in order, separated by commas.
# We use a special character for infinity for clear representation.
print(f"{heesch_number_top}, {heesch_number_middle}, {heesch_number_bottom}")