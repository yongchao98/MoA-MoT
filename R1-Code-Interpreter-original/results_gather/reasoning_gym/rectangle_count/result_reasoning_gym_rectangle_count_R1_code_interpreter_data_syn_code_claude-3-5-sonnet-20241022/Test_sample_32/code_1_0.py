# Initialize variables
rectangle_count = 0
overlap_count = 0
prev_line = ""
in_rectangle = False

# Process input line by line
while True:
    try:
        line = input()
        
        # Count overlaps in this line
        overlap_count += line.count('█')
        
        # Check for rectangle corners
        if '#' in line or '█' in line:
            for i in range(len(line)-1):
                if line[i] in '#█':
                    # Check if it's a top-left corner
                    if (i == 0 or line[i-1] not in '#█') and \
                       (not prev_line or prev_line[i] not in '#█'):
                        rectangle_count += 1
        
        prev_line = line
    except EOFError:
        break

# Calculate final count
# Each overlap point indicates an additional rectangle
total = rectangle_count + overlap_count - 1

print(f"<<<{total}>>>")