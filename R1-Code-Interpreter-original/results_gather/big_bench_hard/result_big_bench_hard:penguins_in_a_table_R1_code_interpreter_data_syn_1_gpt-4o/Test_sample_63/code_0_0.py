# List of penguins with their ages
penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Bernard", "age": 5},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8}
]

# Count penguins with age less than 10
count = sum(1 for penguin in penguins if penguin["age"] < 10)

# Output the count
print(count)