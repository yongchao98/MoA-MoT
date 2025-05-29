def check_statements():
    # Assume Daniel is a sage (True) or a fool (False)
    for daniel_is_sage in [True, False]:
        # Daniel's statement: "Sebastian is a sage and Daniel is a sage"
        daniel_statement = (daniel_is_sage and daniel_is_sage)
        
        # Assume Sebastian is a sage (True) or a fool (False)
        for sebastian_is_sage in [True, False]:
            # Sebastian's statement: "Daniel is a fool if and only if Daniel is a sage"
            sebastian_statement = (not daniel_is_sage) == daniel_is_sage
            
            # Check if both statements are consistent with their roles
            if (daniel_is_sage == daniel_statement) and (sebastian_is_sage == sebastian_statement):
                return daniel_is_sage, sebastian_is_sage

result = check_statements()
print(result)