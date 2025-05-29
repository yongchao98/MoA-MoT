def check_statements():
    # Possible states: (Matthew, David)
    # True represents saint, False represents sinner
    possibilities = [(True, True), (True, False), (False, True), (False, False)]
    
    for matthew, david in possibilities:
        # Matthew's statement: "Matthew is a saint if and only if David is a saint"
        matthew_statement = (matthew == david)
        
        # David's statement: "If Matthew is a saint then Matthew is a sinner"
        david_statement = (not matthew or not matthew)
        
        # Check if the statements align with their nature
        if matthew == matthew_statement and david == david_statement:
            return matthew, david

result = check_statements()
print(result)