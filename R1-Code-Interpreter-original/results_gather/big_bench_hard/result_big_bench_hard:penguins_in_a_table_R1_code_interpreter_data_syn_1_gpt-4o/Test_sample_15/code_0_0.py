penguins = [
    {"name": "Louis", "age": 7, "height": 50, "weight": 11},
    {"name": "Vincent", "age": 9, "height": 60, "weight": 11},
    {"name": "Gwen", "age": 8, "height": 70, "weight": 15},
    {"name": "James", "age": 12, "height": 90, "weight": 12}
]

count = sum(1 for penguin in penguins if penguin["age"] > 5 and penguin["weight"] > 12)
print(count)