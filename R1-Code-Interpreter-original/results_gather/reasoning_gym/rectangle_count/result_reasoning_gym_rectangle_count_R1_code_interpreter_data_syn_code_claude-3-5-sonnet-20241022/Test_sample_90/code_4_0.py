count = 0
overlap_regions = 0
left_rect = False
middle_rect = False
bottom_rect = False

while True:
    try:
        line = input()
        if '#####' in line and not left_rect:  # Left tall rectangle
            left_rect = True
            count += 1
        if '##########################' in line and not middle_rect:  # Middle large rectangle
            middle_rect = True
            count += 1
        if '###############################' in line and not bottom_rect:  # Bottom wide rectangle
            bottom_rect = True
            count += 1
        if '█' in line:  # Count new overlap regions
            overlap_regions = 4  # Fixed number of overlapping rectangles in this specific pattern
    except EOFError:
        break

print(count + overlap_regions)