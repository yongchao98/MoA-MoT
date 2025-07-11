# Define the lengths of the segments of the white path based on the fractal generation rule.

# The first segment is the original main stem.
segment_1 = 40

# The second segment is the stem of the structure from the first iteration.
# Its length is equal to the length of the branch it replaced.
segment_2 = 20

# The third segment is the stem of the structure from the second iteration.
# To find its length, we first need the length of the branch it replaced.
# The branches from the first iteration are scaled by (20/40) = 0.5.
# So, their length is 20 * 0.5 = 10.
# The stem of the second iteration structure is equal to this branch length.
segment_3 = 10

# The last two segments are the two small branches from the second iteration.
# They are scaled by a factor of (10/40) = 0.25 relative to the original branches.
# Their length is 20 * 0.25 = 5 each.
segment_4 = 5
segment_5 = 5

# Calculate the total length by summing all the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4 + segment_5

print("The total length of the white branches is the sum of its segments:")
print(f"{segment_1} + {segment_2} + {segment_3} + {segment_4} + {segment_5} = {total_length}")