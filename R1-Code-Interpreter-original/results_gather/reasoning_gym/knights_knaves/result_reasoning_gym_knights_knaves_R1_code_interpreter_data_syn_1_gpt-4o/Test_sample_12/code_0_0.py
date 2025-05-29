def check_statements():
    # Assume Samuel is a saint
    samuel_is_saint = True
    jacob_is_sinner = not samuel_is_saint  # Based on Samuel's statement

    # Check Jacob's statement
    if samuel_is_saint:
        jacob_statement = samuel_is_saint  # Jacob says "If Samuel is a saint, then Samuel is a sinner"
    else:
        jacob_statement = not samuel_is_saint

    # If Jacob is a sinner, his statement should be false
    jacob_is_sinner = not jacob_statement

    return samuel_is_saint, jacob_is_sinner

samuel_is_saint, jacob_is_sinner = check_statements()
print(f"Samuel is a {'saint' if samuel_is_saint else 'sinner'}, and Jacob is a {'sinner' if jacob_is_sinner else 'saint'}")