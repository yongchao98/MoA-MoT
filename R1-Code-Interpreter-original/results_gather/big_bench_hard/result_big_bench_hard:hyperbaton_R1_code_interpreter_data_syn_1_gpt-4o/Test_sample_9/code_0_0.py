def check_adjective_order(option_a, option_b):
    # Define the correct order of adjectives
    correct_order = ["size", "shape", "color"]
    
    # Define the adjectives and their types for each option
    adjectives_a = [("enormous", "size"), ("rectangular", "shape"), ("blue", "color")]
    adjectives_b = [("enormous", "size"), ("blue", "color"), ("rectangular", "shape")]
    
    # Check the order of adjectives in each option
    order_a = [adj_type for _, adj_type in adjectives_a]
    order_b = [adj_type for _, adj_type in adjectives_b]
    
    # Compare with the correct order
    is_a_correct = order_a == correct_order
    is_b_correct = order_b == correct_order
    
    return is_a_correct, is_b_correct

# Check which option is correct
option_a_correct, option_b_correct = check_adjective_order(
    "enormous rectangular blue dog", 
    "enormous blue rectangular dog"
)

print("Option A is correct:", option_a_correct)
print("Option B is correct:", option_b_correct)