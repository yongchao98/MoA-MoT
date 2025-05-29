# Define the adjectives and their categories
adjectives = {
    "American": "proper",
    "triangular": "shape"
}

# Define the standard order of adjectives
order = ["quantity", "quality", "size", "age", "shape", "color", "proper", "purpose"]

# Function to determine the order of adjectives
def check_adjective_order(adj1, adj2):
    index1 = order.index(adjectives[adj1])
    index2 = order.index(adjectives[adj2])
    return index1 < index2

# Check the order for the given options
option_A = check_adjective_order("American", "triangular")
option_B = check_adjective_order("triangular", "American")

print("Option A is correct:", option_A)
print("Option B is correct:", option_B)