def check_statements():
    # Daniel's statement: "If Ethan is an angel, then Ethan is a devil."
    # Ethan's statement: "Daniel is an angel if and only if Daniel is a devil."
    
    # Assume Daniel is a devil
    daniel_is_angel = False
    # Assume Ethan is a devil
    ethan_is_angel = False
    
    # Check Daniel's statement
    daniel_statement = not ethan_is_angel or ethan_is_angel  # This should be false if Daniel is a devil
    
    # Check Ethan's statement
    ethan_statement = (daniel_is_angel and not daniel_is_angel) or (not daniel_is_angel and daniel_is_angel)  # This should be false if Ethan is a devil
    
    return daniel_statement, ethan_statement

daniel_statement, ethan_statement = check_statements()
print("Daniel's statement is:", daniel_statement)
print("Ethan's statement is:", ethan_statement)