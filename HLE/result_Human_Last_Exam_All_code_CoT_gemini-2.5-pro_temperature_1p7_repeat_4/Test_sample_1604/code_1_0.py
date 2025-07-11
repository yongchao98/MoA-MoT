import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler from a predefined list.
    """
    paintings = {
        "The Roll Call": {"battle": "Battle of Inkerman", "year": 1854},
        "The 28th Regiment at Quatre Bras": {"battle": "Battle of Quatre Bras", "year": 1815},
        "Balaclava": {"battle": "Battle of Balaclava", "year": 1854},
        "The Defence of Rorke's Drift": {"battle": "Battle of Rorke's Drift", "year": 1879},
        "Scotland for Ever!": {"battle": "Battle of Waterloo", "year": 1815},
        "Floreat Etona!": {"battle": "Battle of Laing's Nek", "year": 1881}
    }

    latest_year = 0
    latest_painting_info = {}

    # Find the painting with the latest year
    for title, info in paintings.items():
        if info["year"] > latest_year:
            latest_year = info["year"]
            latest_painting_info = {
                "title": title,
                "battle": info["battle"],
                "year": info["year"]
            }
            
    # Extract the year and its individual digits for the final print statement
    year_str = str(latest_painting_info['year'])
    y1, y2, y3, y4 = year_str[0], year_str[1], year_str[2], year_str[3]

    print(f"The latest historical battle depicted in a painting by Lady Butler is the '{latest_painting_info['battle']}'.")
    print(f"This battle is shown in the painting '{latest_painting_info['title']}' and occurred in the year {y1}{y2}{y3}{y4}.")


find_latest_battle_painting()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)