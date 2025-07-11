def solve_dialect_riddle():
    """
    Solves a riddle by translating a number from a Cumbrian dialect
    to a Derbyshire dialect.
    """
    
    # These are traditional sheep-counting numbers.
    # The spelling can vary, so we use a standard reference.
    cumbrian_dialect = {
        1: 'Yan',
        2: 'Tyan',
        9: 'Dovera', # Phonetically sounds like 'daoves'
        15: 'Bumfit',
        17: 'Tyan-a-bumfit' # 'tyaan'eboon' is a phonetic version of this (2 + 15)
    }

    derbyshire_dialect = {
        1: 'Yan',
        2: 'Tan',
        15: 'Bumper',
        17: 'Tan-bumper' # (2 + 15)
    }

    # The original number from the Kirkby Lonsdale (Cumbrian) person.
    cumbrian_word = "tyaan'eboon"
    original_number = 17
    
    # The word for that number in the Derbyshire dialect.
    derbyshire_equivalent_word = derbyshire_dialect[original_number]

    print(f"The term from Kirkby Lonsdale (Cumbrian dialect) is '{cumbrian_word}'.")
    print(f"This term corresponds to the number 17 in the Cumbrian sheep-counting system ({cumbrian_dialect[2]} + {cumbrian_dialect[15]}).")
    print(f"The equation is: 2 ('Tyan') + 15 ('Bumfit') = 17 ('Tyan-a-bumfit').")
    print("\nIn the Derbyshire dialect:")
    print(f"The word for 17 is '{derbyshire_equivalent_word}'.")
    print(f"This is formed from: 2 ('Tan') + 15 ('Bumper') = 17 ('Tan-bumper').")
    print(f"\nTherefore, a Derbyshireman would have said he had '{derbyshire_equivalent_word}'.")

solve_dialect_riddle()