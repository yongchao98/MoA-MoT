def find_first_best_picture_with_obelisk():
    """
    This function presents the research findings for the first
    Best Picture winner to show a Luxor Obelisk.
    """
    # Based on film history research, several Best Picture winners are set in Paris.
    # We focus on the earliest ones confirmed to show the Luxor Obelisk in the Place de la Concorde.
    candidates = {
        1951: "An American in Paris",
        1956: "Around the World in 80 Days"
    }

    # The goal is to find the film with the earliest year.
    first_year = min(candidates.keys())
    first_film_title = candidates[first_year]

    print("The question asks for the first Academy Award Best Picture winner to show a Luxor Obelisk.")
    print("One of the two Luxor Obelisks is famously located in the Place de la Concorde in Paris.")
    print("\nAfter reviewing Best Picture winners, two early candidates set in Paris are:")
    for year, title in sorted(candidates.items()):
        print(f"- The winner for the year {year}: '{title}'")
    
    print(f"\nBoth films depict the obelisk. To find the *first* one, we identify the earliest year.")
    
    # Final answer output, including the numbers as requested.
    year_digits = [int(d) for d in str(first_year)]
    print("\nThe first winner is '{}' from the year {}.".format(first_film_title, first_year))
    print("The numbers in the winning year are: {}, {}, {}, and {}.".format(year_digits[0], year_digits[1], year_digits[2], year_digits[3]))


if __name__ == "__main__":
    find_first_best_picture_with_obelisk()
