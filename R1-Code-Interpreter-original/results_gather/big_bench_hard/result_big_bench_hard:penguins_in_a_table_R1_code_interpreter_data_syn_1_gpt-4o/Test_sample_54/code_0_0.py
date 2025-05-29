# Initial list of penguins
penguins = [
    {"name": "Louis", "age": 7, "height": 50, "weight": 11},
    {"name": "Bernard", "age": 5, "height": 80, "weight": 13},
    {"name": "Vincent", "age": 9, "height": 60, "weight": 11},
    {"name": "Gwen", "age": 8, "height": 70, "weight": 15}
]

# Add the new penguin
penguins.append({"name": "James", "age": 12, "height": 90, "weight": 12})

# Remove the penguin named Bernard
penguins = [penguin for penguin in penguins if penguin["name"] != "Bernard"]

# Get the name of the last penguin
last_penguin_name = penguins[-1]["name"]

print(last_penguin_name)