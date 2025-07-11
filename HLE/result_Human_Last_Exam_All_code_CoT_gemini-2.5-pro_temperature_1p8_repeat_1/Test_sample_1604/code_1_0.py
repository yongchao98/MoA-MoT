def find_latest_battle():
    """
    Finds the latest historical battle/event depicted in a painting by Lady Butler
    from a predefined list.
    """
    # A list of dictionaries, each containing the painting, the event, and the year of the event.
    paintings_data = [
        {'painting': 'The 28th Regiment at Quatre Bras', 'event': 'Battle of Quatre Bras', 'year': 1815},
        {'painting': 'Scotland Forever!', 'event': 'Battle of Waterloo', 'year': 1815},
        {'painting': 'The Remnants of an Army', 'event': 'Retreat from Kabul', 'year': 1842},
        {'painting': 'The Roll Call', 'event': 'Battle of Inkerman', 'year': 1854},
        {'painting': 'Balaclava', 'event': 'Battle of Balaclava', 'year': 1854},
        {'painting': 'The Defence of Rorke\'s Drift', 'event': 'Battle of Rorke\'s Drift', 'year': 1879},
        {'painting': 'Floreat Etona!', 'event': 'Battle of Laing\'s Nek', 'year': 1881},
        {'painting': 'A Glimpse of the Great Retreat', 'event': 'The Great Retreat', 'year': 1914}
    ]

    latest_event = None
    max_year = 0

    print("Comparing the years of events depicted in Lady Butler's paintings:")

    # Find the maximum year by iterating through the list
    for item in paintings_data:
        event_name = item['event']
        event_year = item['year']
        print(f"- {event_name}: {event_year}")

        if event_year > max_year:
            max_year = event_year
            latest_event = item

    if latest_event:
        print("\n------------------------------------------------------")
        print("To find the latest battle, we compare the years: 1815, 1815, 1842, 1854, 1854, 1879, 1881, 1914.")
        print(f"The latest year is {latest_event['year']}.")
        print(f"The latest historical event depicted is the '{latest_event['event']}', shown in her painting '{latest_event['painting']}'.")
        print("------------------------------------------------------")


find_latest_battle()