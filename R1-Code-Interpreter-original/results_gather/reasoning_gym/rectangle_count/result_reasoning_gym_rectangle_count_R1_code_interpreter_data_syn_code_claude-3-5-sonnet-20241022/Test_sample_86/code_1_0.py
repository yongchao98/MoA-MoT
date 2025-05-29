def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    # Find all continuous segments of '#' or '█' characters
    def find_continuous_segments(row):
        segments = []
        start = None
        for i in range(width):
            if row[i] in '#█':
                if start is None:
                    start = i
            elif start is not None:
                segments.append((start, i - 1))
                start = None
        if start is not None:
            segments.append((start, width - 1))
        return segments

    # Find rectangles starting from each row
    rectangles = 0
    processed = set()  # To avoid counting the same rectangle multiple times
    
    for top in range(height):
        segments = find_continuous_segments(grid[top])
        
        for seg_start, seg_end in segments:
            # Look for matching bottom row
            for bottom in range(top, height):
                # Check if all columns between top and bottom form continuous lines
                valid = True
                for x in range(seg_start, seg_end + 1):
                    for y in range(top, bottom + 1):
                        if grid[y][x] not in '#█':
                            valid = False
                            break
                    if not valid:
                        break
                
                # Check if bottom row matches top row pattern
                if valid:
                    bottom_seg = grid[bottom][seg_start:seg_end + 1]
                    if all(c in '#█' for c in bottom_seg):
                        # Create rectangle signature to avoid duplicates
                        rect_sig = (top, bottom, seg_start, seg_end)
                        if rect_sig not in processed:
                            processed.add(rect_sig)
                            rectangles += 1

    return rectangles

# Create the grid (same grid as before)
grid = [
    # ... (same grid content as before)
]

print(find_rectangles(grid))