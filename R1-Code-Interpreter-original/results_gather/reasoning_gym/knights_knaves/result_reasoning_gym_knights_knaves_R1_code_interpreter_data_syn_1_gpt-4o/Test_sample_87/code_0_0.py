def evaluate_statements():
    # Possible states: (Evelyn, Liam)
    # True represents pioneer, False represents laggard
    possibilities = [(True, True), (True, False), (False, True), (False, False)]
    
    for evelyn, liam in possibilities:
        # Evaluate Evelyn's statement: If Evelyn is a pioneer, then Liam is a laggard
        evelyn_statement = not evelyn or not liam
        
        # Evaluate Liam's statement: Evelyn is a pioneer if and only if Liam is a pioneer
        liam_statement = (evelyn == liam)
        
        # Check if both statements are consistent with their roles
        if (evelyn and evelyn_statement) or (not evelyn and not evelyn_statement):
            if (liam and liam_statement) or (not liam and not liam_statement):
                return (evelyn, liam)

result = evaluate_statements()
print(result)