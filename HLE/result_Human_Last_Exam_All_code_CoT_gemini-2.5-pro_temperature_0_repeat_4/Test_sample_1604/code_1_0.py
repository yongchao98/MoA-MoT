def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler
    from a predefined list.
    """
    # A dictionary of Lady Butler's paintings with their corresponding battle and year.
    # Format: "Painting Title": ("Battle Name", year)
    paintings = {
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Scotland Forever!": ("Battle of Waterloo", 1815),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881),
        "The Camel Corps": ("Battle of Abu Klea", 1885)
    }

    # Find the painting with the latest battle year
    latest_painting_title, (latest_battle_name, latest_year) = max(
        paintings.items(), key=lambda item: item[1][1]
    )

    # Print the result
    print(f"Lady Butler's painting '{latest_painting_title}' depicts the '{latest_battle_name}'.")
    print(f"This battle occurred in the year {latest_year}, making it the latest among her famous works.")
    print(f"The latest battle is the Battle of Abu Klea which happened in 1885.")


if __name__ == "__main__":
    find_latest_battle_painting()