from math import comb

# Total number of rectangles
total_rectangles = comb(3, 2) * comb(2004, 2)

# Rectangles including the shaded column
rectangles_including_shaded = comb(3, 2) * (comb(1002, 2) + comb(1002, 2))

# Rectangles not including the shaded column
rectangles_not_including_shaded = total_rectangles - rectangles_including_shaded

# Probability
probability = rectangles_not_including_shaded / total_rectangles

print(probability)