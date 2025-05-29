def check_statements():
    # Possible scenarios
    scenarios = [
        {"Lily": "hero", "Jack": "villain"},
        {"Lily": "villain", "Jack": "hero"}
    ]
    
    for scenario in scenarios:
        Lily = scenario["Lily"]
        Jack = scenario["Jack"]
        
        # Evaluate Lily's statement
        if Lily == "hero":
            lily_statement = True
        else:
            lily_statement = (Jack == "villain")
        
        # Evaluate Jack's statement
        if Jack == "hero":
            jack_statement = (Lily == "villain")
        else:
            jack_statement = True  # Jack is a villain, so his statement is false
        
        # Check if both statements are consistent
        if lily_statement and jack_statement:
            return scenario

result = check_statements()
print(result)