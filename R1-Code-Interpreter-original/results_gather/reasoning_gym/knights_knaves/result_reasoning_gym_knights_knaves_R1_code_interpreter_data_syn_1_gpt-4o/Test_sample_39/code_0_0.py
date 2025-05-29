def check_statements():
    # Possible states: True for saint, False for sinner
    for isabella in [True, False]:
        for joseph in [True, False]:
            # Isabella's statement: "Joseph is a sinner"
            isabella_statement = not joseph
            
            # Joseph's statement: "Isabella is a saint if and only if Joseph is a saint"
            joseph_statement = (isabella == joseph)
            
            # Check if the statements align with their nature
            if (isabella == isabella_statement) and (joseph == joseph_statement):
                return isabella, joseph

isabella, joseph = check_statements()
print(f"Isabella is a {'saint' if isabella else 'sinner'}, and Joseph is a {'saint' if joseph else 'sinner'}")