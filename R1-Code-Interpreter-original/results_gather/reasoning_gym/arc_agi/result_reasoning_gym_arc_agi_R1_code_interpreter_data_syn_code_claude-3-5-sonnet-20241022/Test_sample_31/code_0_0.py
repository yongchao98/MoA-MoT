def analyze_grid(input_grid):
    # Convert input string to 2D grid
    rows = []
    unique_nums = set()
    for row in input_grid.strip().split('\n'):
        nums = [int(x) for x in row.split()]
        rows.append(nums)
        unique_nums.update(nums)
    
    # Find the background number (most common number)
    background = max(unique_nums, key=lambda x: sum(row.count(x) for row in rows))
    
    # Find non-background numbers and their positions
    numbers = {}
    for i in range(len(rows)):
        for j in range(len(rows[i])):
            if rows[i][j] != background:
                if rows[i][j] not in numbers:
                    numbers[rows[i][j]] = []
                numbers[rows[i][j]].append((i, j))
    
    # Print analysis
    print(f"Background number: {background}")
    print("Non-background numbers and their positions:")
    for num in sorted(numbers.keys()):
        print(f"Number {num}: {numbers[num]}")

# Test input
test_input = """6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 2 2 2 2 2 2 2 2 2 2 6 6
6 6 6 6 2 2 2 2 2 2 2 2 2 2 6 6
6 6 6 6 2 2 2 2 2 2 2 2 2 2 6 6
6 6 6 6 2 2 2 2 2 2 2 2 2 2 6 6
6 6 6 6 2 2 9 9 9 9 2 2 2 2 6 6
6 6 6 6 2 2 9 9 9 9 2 2 2 2 6 6
6 6 6 6 6 6 9 9 9 9 6 6 6 6 6 6
6 6 6 6 6 6 9 9 9 9 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 7 7 7 6 6 6 6 6 6 6 6
6 6 6 6 6 7 7 7 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6 3 3 6 6 6
6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6"""

analyze_grid(test_input)