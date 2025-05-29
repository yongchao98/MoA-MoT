# List of giraffes with their ages
giraffes = [
    {"name": "Jody", "age": 5},
    {"name": "Gladys", "age": 10},
    {"name": "Marian", "age": 2},
    {"name": "Donna", "age": 9}
]

# Count giraffes older than 5 years
count = sum(1 for giraffe in giraffes if giraffe["age"] > 5)

# Output the result
print(count)