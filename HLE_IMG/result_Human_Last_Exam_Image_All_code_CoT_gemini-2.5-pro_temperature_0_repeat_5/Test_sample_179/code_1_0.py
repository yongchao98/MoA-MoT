def find_associated_number():
    """
    This function identifies the piece and calculates the number associated with it from the options.
    """
    composition_year = 1810
    beethoven_birth_year = 1770

    # The manuscript for "Für Elise" is dated April 27, 1810.
    # Beethoven was born in December 1770.
    # In April 1810, he had not yet had his 40th birthday.
    age = composition_year - beethoven_birth_year - 1

    print("The piano piece is Beethoven's 'Für Elise'.")
    print("It was composed in the year 1810.")
    print("Beethoven was born in the year 1770.")
    print("To find his age at the time of composition, we can perform the following calculation:")
    print(f"Age ≈ {composition_year} - {beethoven_birth_year}")
    print("Since the piece was written in April and his birthday was in December, he was still 39.")
    print(f"The associated number is: {age}")

find_associated_number()