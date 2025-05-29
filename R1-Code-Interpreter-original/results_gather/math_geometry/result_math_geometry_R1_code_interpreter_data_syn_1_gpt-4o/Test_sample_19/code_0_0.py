# Number of vertices in a cube
vertices = 8

# Number of edges in a cube
edges = 12

# Number of face diagonals in a cube
face_diagonals = 6 * 2  # 6 faces, each with 2 diagonals

# Number of space diagonals in a cube
space_diagonals = (vertices * (vertices - 1) // 2) - edges - face_diagonals

# Total number of diagonals
total_diagonals = face_diagonals + space_diagonals

print(total_diagonals)