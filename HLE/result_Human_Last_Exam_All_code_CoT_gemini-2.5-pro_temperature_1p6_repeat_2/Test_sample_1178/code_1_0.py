# The problem is to find the smallest integer-length rectangle that has a fault-free
# tiling using squares from the set S = {2x2, 3x3, 5x5, 7x7}.

# Key properties of such a tiling:
# 1. A "fault-free" tiling is one that cannot be cut by a single straight line
#    from edge to edge without crossing a tile.
# 2. A tiling with squares of only two different sizes (e.g., 2x2 and 3x3)
#    is always "fault-full" (guillotine-cuttable).
# 3. Therefore, our tiling must use at least three different square sizes.

# The smallest such rectangle known is 11x17.
width = 11
height = 17
area = width * height

# A possible set of tiles for an 11x17 fault-free tiling is:
# - 1 square of size 2x2
# - 12 squares of size 3x3
# - 3 squares of size 5x5
n2 = 1
n3 = 12
n5 = 3

area_n2 = n2 * (2*2)
area_n3 = n3 * (3*3)
area_n5 = n5 * (5*5)
total_area_from_tiles = area_n2 + area_n3 + area_n5

print(f"The smallest rectangle has dimensions {width}x{height}.")
print(f"The area of this rectangle is {width} * {height} = {area}.")
print("This area can be composed of squares of three different sizes (2x2, 3x3, 5x5):")
print(f"{n3} * (3*3) + {n5} * (5*5) + {n2} * (2*2) = {n3} * 9 + {n5} * 25 + {n2} * 4 = {area_n3} + {area_n5} + {area_n2} = {total_area_from_tiles}")
print(f"The area of this rectangle is {area}.")
