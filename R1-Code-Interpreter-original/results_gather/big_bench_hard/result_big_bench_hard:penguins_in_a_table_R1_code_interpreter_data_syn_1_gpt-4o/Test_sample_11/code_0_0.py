# List of giraffes with their attributes
giraffes = [
    {"name": "Jody", "age": 5, "height": 430, "weight": 620},
    {"name": "Gladys", "age": 10, "height": 420, "weight": 590},
    {"name": "Marian", "age": 2, "height": 310, "weight": 410},
    {"name": "Donna", "age": 9, "height": 440, "weight": 650}
]

# Count giraffes older than 5 years
count = sum(1 for giraffe in giraffes if giraffe["age"] > 5)

# Output the result
print(count)