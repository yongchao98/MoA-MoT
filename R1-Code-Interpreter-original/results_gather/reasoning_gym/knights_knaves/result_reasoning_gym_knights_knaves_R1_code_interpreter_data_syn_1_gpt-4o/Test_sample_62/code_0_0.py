def check_statements():
    # Zoey's statement: "If Joseph is a sage, then Joseph is a fool."
    # Joseph's statement: "Zoey is a fool and Zoey is a sage."
    
    # Assume Zoey is a sage
    zoey_is_sage = True
    # If Zoey is a sage, her statement must be true, so Joseph must be a fool
    joseph_is_fool = True
    
    # Check Joseph's statement
    # Joseph's statement is false if he is a fool
    joseph_statement_is_false = joseph_is_fool
    
    # If both conditions hold, our assumption is correct
    if zoey_is_sage and joseph_statement_is_false:
        return "Zoey is a sage, and Joseph is a fool"
    else:
        return "The assumption is incorrect"

print(check_statements())