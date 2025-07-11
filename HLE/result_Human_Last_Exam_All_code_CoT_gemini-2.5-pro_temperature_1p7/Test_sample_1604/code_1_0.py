import operator

def find_latest_battle():
    """
    Identifies the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # A dictionary mapping paintings to the battle depicted and its year.
    battles = {
        "'The Remnants of an Army'": {"event": "1842 retreat from Kabul", "year": 1842},
        "'The Roll Call'": {"event": "Crimean War (specifically, after the Battle of Inkerman)", "year": 1854},
        "'The Defence of Rorke's Drift'": {"event": "Battle of Rorke's Drift", "year": 1879},
        "'Scotland for Ever!'": {"event": "Battle of Waterloo", "year": 1815},
        "'Floreat Etona!'": {"event": "First Boer War (inspired by the Battle of Laing's Nek)", "year": 1881},
        "'After the Battle'": {"event": "Battle of Tel El Kebir", "year": 1882},
        "'The Camel Corps'": {"event": "Nile Expedition (part of the Mahdist War)", "year": 1885}
    }

    # Find the painting with the latest year.
    # The max function's key argument specifies a function to be called on each list element prior to making comparisons.
    # In this case, we're comparing the 'year' value of each dictionary entry.
    latest_painting_name, latest_info = max(battles.items(), key=lambda item: item[1]['year'])

    print(f"To find the latest battle depicted by Lady Butler, we compare the years of the events shown in her paintings:")
    for painting, info in battles.items():
        print(f"- {painting}: {info['event']} in {info['year']}")

    print("\nComparing the years...")
    years = [info['year'] for info in battles.values()]
    print(f"{' > '.join(map(str, sorted(years, reverse=True)))}")

    print(f"\nThe latest year is {latest_info['year']}.")
    print(f"This corresponds to the event '{latest_info['event']}', depicted in the painting {latest_painting_name}.")

find_latest_battle()