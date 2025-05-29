def check_adjective_order(option):
    order = {
        "size": ["midsize"],
        "shape": ["prismlike"],
        "color": ["orange"],
        "material": ["iron"],
        "purpose": ["snorkeling"]
    }
    
    adjectives = option.split()
    order_indices = []
    
    for adj in adjectives:
        for category, words in order.items():
            if adj in words:
                order_indices.append(list(order.keys()).index(category))
    
    return order_indices == sorted(order_indices)

option_a = "snorkeling iron prismlike orange midsize knife"
option_b = "midsize prismlike orange iron snorkeling knife"

result_a = check_adjective_order(option_a)
result_b = check_adjective_order(option_b)

print("Option A is correct:", result_a)
print("Option B is correct:", result_b)