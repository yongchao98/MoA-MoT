# List of giraffes with their ages
giraffes = [
    ("Jody", 5),
    ("Gladys", 10),
    ("Marian", 2),
    ("Donna", 9)
]

# Find the giraffe with the maximum age
oldest_giraffe = max(giraffes, key=lambda giraffe: giraffe[1])

# Output the name of the oldest giraffe
print(oldest_giraffe[0])