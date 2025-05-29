def is_pioneer(statement, is_pioneer):
    # If the person is a pioneer, the statement must be true
    return statement if is_pioneer else not statement

# Chloe's statement: If David is a pioneer, then David is a laggard
# This is equivalent to: not (David is a pioneer and David is a pioneer)
chloe_statement = lambda david_is_pioneer: not (david_is_pioneer and david_is_pioneer)

# David's statement: Chloe is a pioneer and David is a pioneer
david_statement = lambda chloe_is_pioneer, david_is_pioneer: chloe_is_pioneer and david_is_pioneer

# Check all combinations
results = []
for chloe_is_pioneer in [True, False]:
    for david_is_pioneer in [True, False]:
        chloe_truth = is_pioneer(chloe_statement(david_is_pioneer), chloe_is_pioneer)
        david_truth = is_pioneer(david_statement(chloe_is_pioneer, david_is_pioneer), david_is_pioneer)
        if chloe_truth and david_truth:
            results.append((chloe_is_pioneer, david_is_pioneer))

print(results)