def find_latest_butler_battle():
    """
    Identifies the latest historical battle depicted in a well-known painting by Lady Butler.
    """
    # A dictionary mapping the subject of Lady Butler's paintings to the year of the historical event.
    battles_depicted = {
        "The Battle of Quatre Bras (from 'The 28th Regiment at Quatre Bras')": 1815,
        "The Battle of Waterloo (from 'Scotland Forever!')": 1815,
        "The Battle of Inkerman (from 'The Roll Call')": 1854,
        "The Battle of Rorke's Drift (from 'The Defence of Rorke's Drift')": 1879,
        "The Battle of Laing's Nek (from 'Floreat Etona!')": 1881,
        "The Battle of Tel El Kebir (from 'Tel El Kebir')": 1882,
        "The Gordon Relief Expedition (from 'The Camel Corps')": 1885
    }

    # Find the event with the latest year
    latest_event_name = max(battles_depicted, key=battles_depicted.get)
    latest_year = battles_depicted[latest_event_name]

    print(f"The latest historical event from a selection of Lady Butler's famous paintings is:")
    print(f"{latest_event_name}, which occurred in the year {latest_year}.")

find_latest_butler_battle()