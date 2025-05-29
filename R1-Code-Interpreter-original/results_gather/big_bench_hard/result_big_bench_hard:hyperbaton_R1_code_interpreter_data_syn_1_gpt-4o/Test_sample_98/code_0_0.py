# Define the adjectives and their categories
adjectives = {
    "big": "size",
    "square": "shape"
}

# Define the correct order of adjectives
correct_order = ["size", "shape"]

# Check the order of the given options
option_A = ["big", "square"]
option_B = ["square", "big"]

# Function to check if the order is correct
def check_order(option, adjectives, correct_order):
    order = [adjectives[adj] for adj in option]
    return order == correct_order

# Check both options
result_A = check_order(option_A, adjectives, correct_order)
result_B = check_order(option_B, adjectives, correct_order)

print("Option A is correct:", result_A)
print("Option B is correct:", result_B)