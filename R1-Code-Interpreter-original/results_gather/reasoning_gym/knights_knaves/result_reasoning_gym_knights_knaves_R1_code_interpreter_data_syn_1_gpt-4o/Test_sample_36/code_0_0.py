def determine_roles():
    # Assume Abigail is a sage
    abigail_is_sage = True
    david_is_sage = False  # Based on Abigail's statement, this should be false if Abigail is a sage

    # Check Abigail's statement: "David is a sage if and only if David is a fool"
    if abigail_is_sage:
        # This statement is a contradiction, so Abigail cannot be a sage
        abigail_is_sage = False
        david_is_sage = True  # If Abigail is a fool, David's statement must be true

    # Check David's statement: "Abigail is a fool"
    if david_is_sage:
        abigail_is_fool = True
    else:
        abigail_is_fool = False

    return abigail_is_sage, abigail_is_fool, david_is_sage

abigail_is_sage, abigail_is_fool, david_is_sage = determine_roles()
print(f"Abigail is a {'sage' if abigail_is_sage else 'fool'}, and David is a {'sage' if david_is_sage else 'fool'}")