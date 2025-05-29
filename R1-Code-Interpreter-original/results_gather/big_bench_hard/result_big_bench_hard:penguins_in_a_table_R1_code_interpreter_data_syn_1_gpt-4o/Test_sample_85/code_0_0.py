# Initial list of penguins
penguins = [
    {"name": "Louis", "age": 7, "height": 50, "weight": 11},
    {"name": "Bernard", "age": 5, "height": 80, "weight": 13},
    {"name": "Vincent", "age": 9, "height": 60, "weight": 11},
    {"name": "Gwen", "age": 8, "height": 70, "weight": 15}
]

# Add James
penguins.append({"name": "James", "age": 12, "height": 90, "weight": 12})

# Remove Bernard
penguins = [penguin for penguin in penguins if penguin["name"] != "Bernard"]

# Find the penguin with height 60 cm
penguin_60_cm = next((penguin["name"] for penguin in penguins if penguin["height"] == 60), None)

print(penguin_60_cm)