def find_kanji_radical():
    """
    This function identifies the radical (bushu) of a specific Kanji character
    and prints its name in Hepburn romanization.
    """
    # The character we are analyzing.
    target_character = '謄'

    # A simple knowledge base mapping a few kanji to their radicals' names.
    # The radical for 謄 is 月, which in this context is a variant of the
    # 'meat' radical (肉) called 'Nikuzuki' (にくづき).
    radical_database = {
        '謄': 'Nikuzuki',
        '媛': 'Onnahen'  # The example from the user's prompt
    }

    # Retrieve the radical's name from our database.
    radical_name = radical_database.get(target_character, "Radical not found.")

    # Print the result in a clear format.
    print(f"The bushu / radical in the Japanese character '{target_character}' is '{radical_name}'.")

if __name__ == "__main__":
    find_kanji_radical()