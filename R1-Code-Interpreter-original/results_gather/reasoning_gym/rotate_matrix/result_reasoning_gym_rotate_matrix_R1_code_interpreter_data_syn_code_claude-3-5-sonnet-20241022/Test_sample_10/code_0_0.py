def rotate_180(matrix):
    n = len(matrix)
    # Create a flat list of all numbers
    flat = [num for row in matrix for num in row]
    # Reverse the list
    flat.reverse()
    # Reconstruct the matrix
    result = []
    for i in range(0, len(flat), n):
        result.append(flat[i:i+n])
    
    # Format the output string
    output = ""
    for row in result:
        output += " ".join(map(str, row)) + "\n"
    output = output.strip()
    
    print(f"<<<{output}>>>")

# Input matrix
matrix = [
    [0, 4, 3],
    [3, 2, 1],
    [8, 1, 9]
]

rotate_180(matrix)