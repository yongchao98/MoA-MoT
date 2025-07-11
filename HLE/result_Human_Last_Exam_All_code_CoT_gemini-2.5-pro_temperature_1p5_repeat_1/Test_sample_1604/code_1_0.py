import sys
# Redirect stdout to a file to prevent gradio from complaining about printing to stdout
# This is a workaround for a specific environment issue and may not be necessary in all contexts.
# sys.stdout = open('output.txt', 'w')

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # A dictionary of Lady Butler's paintings, the battles they depict, and the year of the battle.
    # Note: For campaigns spanning multiple years, the end year is used for comparison.
    paintings = {
        'Quatre Bras': {'battle': 'Battle of Quatre Bras', 'year': 1815},
        'Scotland for Ever!': {'battle': 'Battle of Waterloo', 'year': 1815},
        'The Remnants of an Army': {'battle': 'Retreat from Kabul', 'year': 1842},
        'The Roll Call': {'battle': 'Battle of Inkerman', 'year': 1854},
        'The Defence of Rorke\'s Drift': {'battle': 'Battle of Rorke\'s Drift', 'year': 1879},
        'Floreat Etona!': {'battle': 'Battle of Laing\'s Nek', 'year': 1881},
        'Tel-el-Kebir': {'battle': 'Battle of Tel El Kebir', 'year': 1882},
        'The Camel Corps': {'battle': 'Gordon Relief Expedition', 'year': 1885},
        'A Lament in the Desert': {'battle': 'Sinai and Palestine campaign (WWI)', 'year': 1918},
    }

    # Find the painting corresponding to the latest battle year
    latest_painting_title = max(paintings, key=lambda p: paintings[p]['year'])
    
    latest_battle_info = paintings[latest_painting_title]
    battle_name = latest_battle_info['battle']
    battle_year = latest_battle_info['year']

    print(f"The latest historical event depicted in a known painting by Lady Butler is from the year {battle_year}.")
    print(f"The painting is '{latest_painting_title}' and it depicts a scene from the {battle_name}.")

# Execute the function
find_latest_battle()

# sys.stdout.close()
# sys.stdout = sys.__stdout__
