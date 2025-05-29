# Define the houses and their attributes
houses = [
    {"name": "Bob", "lunch": "soup", "color": "white", "phone": "OnePlus 9"},
    {"name": "Arnold", "lunch": None, "color": "purple", "phone": "Samsung Galaxy S21"},
    {"name": "Alice", "lunch": "pizza", "color": "blue", "phone": "Huawei P50"},
    {"name": "Carol", "lunch": "grilled cheese", "color": "brown", "phone": None}
]

# Check the conditions
conditions_met = (
    houses[0]["lunch"] == "soup" and
    houses[0]["color"] == "white" and
    houses[0]["phone"] == "OnePlus 9" and
    houses[1]["name"] == "Arnold" and
    houses[1]["color"] == "purple" and
    houses[1]["phone"] == "Samsung Galaxy S21" and
    houses[2]["name"] == "Alice" and
    houses[2]["color"] == "blue" and
    houses[2]["phone"] == "Huawei P50" and
    houses[3]["name"] == "Carol" and
    houses[3]["color"] == "brown" and
    houses[3]["lunch"] == "grilled cheese"
)

print("All conditions met:", conditions_met)