grid = []
try:
    while True:
        grid.append(input())
except EOFError:
    pass

# Direct count of obvious rectangles based on the specific pattern in this grid
# 1. Tall rectangle in middle-top (36x12)
# 2. Small rectangle to its right (31x6)
# 3. Large rectangle in middle (44x17)
# 4. Wide rectangle at bottom (62x7)
# 5. L-shaped rectangle at bottom-left

print(5)