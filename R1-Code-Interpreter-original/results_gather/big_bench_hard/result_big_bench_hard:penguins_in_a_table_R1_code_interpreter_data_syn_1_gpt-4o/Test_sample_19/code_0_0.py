penguins = [
    {"name": "Louis", "age": 7, "height": 50, "weight": 11},
    {"name": "Vincent", "age": 9, "height": 60, "weight": 11},
    {"name": "Gwen", "age": 8, "height": 70, "weight": 15}
]

youngest_penguin = min(penguins, key=lambda x: x["age"])
print(youngest_penguin["name"])