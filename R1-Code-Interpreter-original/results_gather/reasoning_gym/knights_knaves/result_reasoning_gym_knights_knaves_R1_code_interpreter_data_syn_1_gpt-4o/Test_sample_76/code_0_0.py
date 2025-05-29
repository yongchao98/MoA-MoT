def analyze_statements():
    # Possible states: True for hero, False for villain
    for benjamin in [True, False]:
        for scarlett in [True, False]:
            # Benjamin's statement: If Benjamin is a hero, then Scarlett is a hero
            benjamin_statement = not benjamin or scarlett
            
            # Scarlett's statement: Scarlett is a hero or Benjamin is a hero
            scarlett_statement = scarlett or benjamin
            
            # Check if both statements are consistent with their nature
            if (benjamin == benjamin_statement) and (scarlett == scarlett_statement):
                return benjamin, scarlett

benjamin_is_hero, scarlett_is_hero = analyze_statements()
print(f"Benjamin is a {'hero' if benjamin_is_hero else 'villain'}, and Scarlett is a {'hero' if scarlett_is_hero else 'villain'}")