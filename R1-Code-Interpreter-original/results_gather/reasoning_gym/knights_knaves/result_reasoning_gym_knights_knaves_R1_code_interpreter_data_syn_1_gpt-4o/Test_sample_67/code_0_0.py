def check_statements():
    # Possible scenarios
    scenarios = [
        {"Matthew": "pioneer", "Ava": "laggard"},
        {"Matthew": "laggard", "Ava": "pioneer"}
    ]
    
    for scenario in scenarios:
        matthew_is_pioneer = scenario["Matthew"] == "pioneer"
        ava_is_laggard = scenario["Ava"] == "laggard"
        
        # Matthew's statement: "Matthew is a pioneer if and only if Ava is a laggard"
        matthew_statement = matthew_is_pioneer == ava_is_laggard
        
        # Ava's statement: "If Matthew is a laggard then Matthew is a pioneer"
        ava_statement = not matthew_is_pioneer or matthew_is_pioneer
        
        if matthew_statement and not ava_statement:
            return scenario

result = check_statements()
print(f"Matthew is a {result['Matthew']}, and Ava is a {result['Ava']}")