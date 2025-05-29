def analyze_brackets(s):
    count_open_square = s.count('[')
    count_close_square = s.count(']')
    count_open_angle = s.count('<')
    count_close_angle = s.count('>')
    
    print(f"Opening square brackets '[': {count_open_square}")
    print(f"Closing square brackets ']': {count_close_square}")
    print(f"Opening angle brackets '<': {count_open_angle}")
    print(f"Closing angle brackets '>': {count_close_angle}")
    
    missing = ']' * (count_open_square - count_close_square)
    print(f"\nMissing brackets: {missing}")
    print(f"Complete sequence: {s + missing}")

input_str = "[ < > ] [ [ < > ]"
analyze_brackets(input_str)