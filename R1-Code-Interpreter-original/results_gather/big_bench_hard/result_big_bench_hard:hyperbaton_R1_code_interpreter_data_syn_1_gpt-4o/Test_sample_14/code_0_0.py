# Define the adjectives and their order
adjective_order = {
    "quantity": 1,
    "quality": 2,
    "size": 3,
    "age": 4,
    "shape": 5,
    "color": 6,
    "proper": 7,
    "purpose": 8
}

# Define the adjectives in the options
option_a = ["nice", "American"]
option_b = ["American", "nice"]

# Define a function to get the order of adjectives
def get_adjective_order(adjectives):
    order = []
    for adj in adjectives:
        if adj == "nice":
            order.append(adjective_order["quality"])
        elif adj == "American":
            order.append(adjective_order["proper"])
    return order

# Check the order for both options
order_a = get_adjective_order(option_a)
order_b = get_adjective_order(option_b)

# Determine which option has the correct order
correct_option = "A" if order_a == sorted(order_a) else "B"

print(correct_option)