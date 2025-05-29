def check_adjective_order(sentence):
    adjectives = sentence.split()[:-1]  # Exclude the noun
    order = ["quantity", "quality", "size", "age", "shape", "color", "proper", "purpose"]
    adjective_types = {"blue": "color", "driving": "purpose"}
    
    indices = [order.index(adjective_types[adj]) for adj in adjectives if adj in adjective_types]
    
    return indices == sorted(indices)

# Check both options
option_a = "driving blue car"
option_b = "blue driving car"

result_a = check_adjective_order(option_a)
result_b = check_adjective_order(option_b)

print("Option A is correct:", result_a)
print("Option B is correct:", result_b)