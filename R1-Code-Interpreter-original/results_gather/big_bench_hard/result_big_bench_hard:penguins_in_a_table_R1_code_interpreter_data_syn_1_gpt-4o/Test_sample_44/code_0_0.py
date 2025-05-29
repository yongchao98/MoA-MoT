# List of penguins with their ages
penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Bernard", "age": 5},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8},
    {"name": "James", "age": 12}
]

# Count penguins older than 8 years
count = sum(1 for penguin in penguins if penguin["age"] > 8)

# Output the result
print(count)