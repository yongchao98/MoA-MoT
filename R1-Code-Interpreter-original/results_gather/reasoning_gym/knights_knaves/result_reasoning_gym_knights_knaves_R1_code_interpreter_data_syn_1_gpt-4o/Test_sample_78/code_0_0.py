def check_statements():
    # Assume Sebastian is a sage (True) or a fool (False)
    # Assume Benjamin is a sage (True) or a fool (False)
    for sebastian in [True, False]:
        for benjamin in [True, False]:
            # Sebastian's statement: "Benjamin is a fool"
            sebastian_statement = not benjamin
            
            # Benjamin's statement: "Sebastian is a sage if and only if Benjamin is a sage"
            benjamin_statement = (sebastian == benjamin)
            
            # Check if the statements align with their nature
            if (sebastian == sebastian_statement) and (benjamin == benjamin_statement):
                return f"Sebastian is a {'sage' if sebastian else 'fool'}, and Benjamin is a {'sage' if benjamin else 'fool'}"

print(check_statements())