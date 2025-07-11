def solve_dialect_riddle():
    """
    Solves a riddle by translating a number between Cumbrian and
    Derbyshire sheep-counting dialects.
    """
    # 1. Explanation of the dialect terms.
    # The terms are from Cumbrian (Kirkby Lonsdale) sheep counting.
    cumbrian_term_past = "tyaan'eboon"
    cumbrian_term_present = "daoves"

    # 2. Decoding the numbers.
    # 'Tyan' is 2, and 'eboon' (aboon) means 'above', usually above 10.
    # So, 'tyaan'eboon' represents the number 12.
    # 'Daoves' is a local variant for 10.
    value_past = 12
    value_present = 10

    # The question asks what the person "used to have had", which is the first number.
    target_value = value_past

    # 3. Finding the Derbyshire equivalent.
    # In the Derbyshire dialect: 'Tan' = 2, 'Dik' = 10.
    # Therefore, 12 is 'Tan-a-dik' (two-and-ten).
    derbyshire_term = "Tan-a-dik"
    derbyshire_value = 12

    # 4. Printing the explanation and the result.
    print("This riddle requires translating a number from a Cumbrian dialect to a Derbyshire dialect.")
    print(f"The term 'tyaan\'eboon' from Kirkby Lonsdale means 'two-above-ten', which is the number {value_past}.")
    print(f"The person now has 'daoves' worth, which is {value_present}.")
    print(f"\nThe number he 'used to have had' is therefore {target_value}.")
    print("\nIn the Derbyshire sheep-counting system, the number 12 is called 'Tan-a-dik'.")
    print("\nTherefore, the final answer is:")
    print(f"Cumbrian: {cumbrian_term_past} ({value_past}) = Derbyshire: {derbyshire_term} ({derbyshire_value})")

solve_dialect_riddle()