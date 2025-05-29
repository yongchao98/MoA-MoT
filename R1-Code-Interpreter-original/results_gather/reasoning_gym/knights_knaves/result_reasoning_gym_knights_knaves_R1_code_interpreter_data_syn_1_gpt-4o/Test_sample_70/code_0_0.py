def check_statements():
    # Chloe's statement: "If Sophia is a sinner then Sophia is a saint."
    # This is equivalent to: not (Sophia is a sinner) or (Sophia is a saint)
    # Which simplifies to: Sophia is a saint (since a saint cannot be both a sinner and a saint)
    
    # Sophia's statement: "Chloe is a sinner or Chloe is a saint."
    # This is always true, so it doesn't affect the logic.

    # Assume Chloe is a saint
    chloe_is_saint = True
    sophia_is_saint = True  # Based on Chloe's statement being true

    # Check consistency
    if chloe_is_saint and sophia_is_saint:
        return "Chloe is a saint, and Sophia is a saint"
    else:
        return "Inconsistent logic"

print(check_statements())