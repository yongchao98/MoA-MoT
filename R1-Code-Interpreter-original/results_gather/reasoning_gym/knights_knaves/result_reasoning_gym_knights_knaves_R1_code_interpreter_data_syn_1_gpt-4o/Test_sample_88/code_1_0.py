def is_knight_or_knave(emily_knight, abigail_knight):
    # Emily's statement: "Abigail is a knight or Emily is a knight"
    emily_statement = abigail_knight or emily_knight
    
    # Abigail's statement: "Emily is a knave if and only if Abigail is a knight"
    abigail_statement = (not emily_knight) == abigail_knight
    
    # Check if the statements are consistent with their roles
    emily_truth = emily_statement == emily_knight
    abigail_truth = abigail_statement == abigail_knight
    
    return emily_truth and abigail_truth

# Check all possible combinations
results = []
for emily_knight in [True, False]:
    for abigail_knight in [True, False]:
        if is_knight_or_knave(emily_knight, abigail_knight):
            results.append((emily_knight, abigail_knight))

print(results)