def check_statements():
    # Ella's statement: "Elizabeth is a sinner if and only if Elizabeth is a saint."
    # This is a contradiction, so Ella must be a sinner.
    ella_is_sinner = True

    # Elizabeth's statement: "Elizabeth is a saint and Ella is a saint."
    # Since Ella is a sinner, this statement is false, so Elizabeth must be a sinner.
    elizabeth_is_sinner = True

    return ella_is_sinner, elizabeth_is_sinner

ella_is_sinner, elizabeth_is_sinner = check_statements()
print(f"Ella is a {'sinner' if ella_is_sinner else 'saint'}, and Elizabeth is a {'sinner' if elizabeth_is_sinner else 'saint'}")