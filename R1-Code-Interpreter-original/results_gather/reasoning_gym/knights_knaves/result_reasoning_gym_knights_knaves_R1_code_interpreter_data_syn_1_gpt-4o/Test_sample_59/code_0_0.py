def check_statements():
    # Possible scenarios
    scenarios = [
        {"Avery": "hero", "Aiden": "hero"},
        {"Avery": "hero", "Aiden": "villain"},
        {"Avery": "villain", "Aiden": "hero"},
        {"Avery": "villain", "Aiden": "villain"}
    ]
    
    for scenario in scenarios:
        Avery = scenario["Avery"]
        Aiden = scenario["Aiden"]
        
        # Avery's statement: "Avery is a hero."
        avery_statement = (Avery == "hero")
        
        # Aiden's statement: "If Aiden is a hero, then Avery is a hero."
        aiden_statement = not (Aiden == "hero" and Avery != "hero")
        
        # Check if both statements are consistent with their roles
        if (Avery == "hero" and avery_statement) or (Avery == "villain" and not avery_statement):
            if (Aiden == "hero" and aiden_statement) or (Aiden == "villain" and not aiden_statement):
                return f"Avery is a {Avery}, and Aiden is a {Aiden}"

print(check_statements())