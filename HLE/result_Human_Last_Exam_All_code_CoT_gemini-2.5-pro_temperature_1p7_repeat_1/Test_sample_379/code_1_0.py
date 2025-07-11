import math

# The side length of the cube is r.
# We can use r=1.0 since it will cancel out in the final ratio.
r = 1.0

# 1. Calculate the length of the locus C on each relevant face.
# On each of the 2 faces sharing the edge with P, the arc length is (r * pi / 3).
len_face_P = r * math.pi / 3
# On each of the 2 adjacent "side" faces, the total arc length is (2 * r * pi / 3).
len_face_adj = 2 * r * math.pi / 3

# 2. Calculate the total length of C by summing the contributions from the 4 faces.
# (The other 2 faces are too far away, contributing 0 length).
total_length_C = 2 * len_face_P + 2 * len_face_adj

# 3. Calculate the divisor as per the problem statement.
divisor = 2 * math.pi * r

# 4. Compute the ratio and express as a whole number percentage.
ratio = total_length_C / divisor
percentage = int(round(ratio * 100))

# 5. Output the numbers in the final equation as requested.
print("The total length of C is the sum of lengths on 4 faces.")
# The "2*(...)" terms represent the two faces of each type
print(f"Total Length = 2 * (r * pi / 3) + 2 * (2 * r * pi / 3)")
print(f"Total Length = (2 * pi * r / 3) + (4 * pi * r / 3) = 2 * pi * r")
print("\nThe final equation is: (Total Length of C) / (2 * pi * r) * 100%")
# Showing the equation with values (where r cancels out)
# (2 * pi) / (2 * pi) * 100
print(f"({total_length_C:.4f}) / ({divisor:.4f}) * 100% = {percentage}%")
<<<100>>>