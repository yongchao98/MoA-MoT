# List of penguins after removing Bernard
penguins = [
    {"name": "Louis", "age": 7, "height": 50, "weight": 11},
    {"name": "Vincent", "age": 9, "height": 60, "weight": 11},
    {"name": "Gwen", "age": 8, "height": 70, "weight": 15}
]

# Count penguins more than 5 years old and weight less than 12 kg
count = sum(1 for penguin in penguins if penguin["age"] > 5 and penguin["weight"] < 12)

print(count)