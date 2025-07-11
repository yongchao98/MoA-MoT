import sys

def find_latest_battle_painting():
    """
    Finds and prints the details of the Lady Butler painting
    that depicts the most recent historical battle from a known list.
    """
    # A list of Lady Butler's battle paintings and the years of the events depicted.
    paintings = [
        {'title': 'Steady the Drums and Fifes', 'battle': 'Battle of Albuera', 'year': 1811},
        {'title': 'Scotland for Ever!', 'battle': 'Battle of Waterloo', 'year': 1815},
        {'title': 'The Roll Call', 'battle': 'Battle of Inkerman', 'year': 1854},
        {'title': 'The Defence of Rorke\'s Drift', 'battle': 'Battle of Rorke\'s Drift', 'year': 1879},
        {'title': 'Floreat Etona!', 'battle': 'Battle of Laing\'s Nek', 'year': 1881},
        {'title': 'The Camel Corps', 'battle': 'Gordon Relief Expedition (Battle of Abu Klea)', 'year': 1885}
    ]

    if not paintings:
        print("The list of paintings is empty.")
        return

    # Assume the first painting is the latest to start
    latest_battle_painting = paintings[0]

    # Find the painting with the highest year
    for painting in paintings:
        if painting['year'] > latest_battle_painting['year']:
            latest_battle_painting = painting

    # Output the result
    title = latest_battle_painting['title']
    battle = latest_battle_painting['battle']
    year = latest_battle_painting['year']

    print(f"The latest historical battle depicted in a painting by Lady Butler is from the year {year}.")
    print(f"The painting is titled '{title}' and it shows events from the {battle}.")

find_latest_battle_painting()