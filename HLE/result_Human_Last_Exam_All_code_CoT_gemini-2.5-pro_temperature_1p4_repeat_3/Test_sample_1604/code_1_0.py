def find_latest_butler_battle():
    """
    This function identifies the latest historical battle depicted in a painting by Lady Butler.
    It uses a predefined dictionary of her famous battle paintings and their corresponding battle dates.
    """
    # A dictionary where keys are battle names and values are tuples containing the year and the painting's title.
    battles_info = {
        "Battle of Albuera": (1811, "'Steady the Drums and Fifes'"),
        "Battle of Waterloo": (1815, "'Scotland for Ever!'"),
        "Battle of Balaclava": (1854, "'Balaclava'"),
        "Battle of Inkerman": (1854, "'The Roll Call'"),
        "Battle of Rorke's Drift": (1879, "'The Defence of Rorke's Drift'"),
        "Battle of Laing's Nek": (1881, "'Floreat Etona!'"),
        "Battle of Tel El Kebir": (1882, "'Tel-el-Kebir'"),
        "Battle of Abu Klea": (1885, "'The Camel Corps'"),
        "Affair of Agagia": (1916, "'The Dorset Yeomanry at Agagia, 26th February 1916'")
    }

    # Find the entry with the maximum year.
    # The `key` for max() is a lambda function that tells it to look at the year (the first element of the tuple value).
    latest_battle_name, (latest_year, painting_title) = max(battles_info.items(), key=lambda item: item[1][0])

    # Print the result in a clear, single sentence.
    print(f"The latest battle depicted by Lady Butler is the '{latest_battle_name}' which occurred in {latest_year}, as shown in her painting {painting_title}.")

find_latest_butler_battle()