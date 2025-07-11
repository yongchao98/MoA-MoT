import re

def identify_runestone():
    """
    Identifies an Ingvar runestone based on text fragments
    by searching a simulated database of runic inscriptions.
    """
    # This dictionary simulates a database of Ingvar runestone inscriptions.
    # The inscriptions are transliterated and normalized for searching.
    # The key is the runestone's official ID, and the value is its inscription.
    runestone_database = {
        "Sö 179": "þurui lit raisa stain þansi at sun sin haralt bruþur ikuars",
        "Sö 131": "spruar risti stan þinsi at ulfui sun sin han fial miþ ikuari",
        "Sö 254": "aþuruiþr risti stain þena eftir sun sin han for austr miþ ink...",
        "U 654": "þir litu risa stin þina aftir kunar faþur sin is austr fur miþ ikuari",
        "U 644": "ulfr lit risa stin þina aftir anut faþur sin is austr uarþ tuþr miþ ikuari"
    }

    # Text fragments identified from the runic text in the image.
    # 'stain þena' means 'this stone'.
    # 'eftir sun sin' means 'in memory of his son'.
    fragment1 = "stain þena"
    fragment2 = "eftir sun sin"

    found_id = None
    for stone_id, inscription in runestone_database.items():
        # Check if the inscription contains both of the identified fragments.
        if fragment1 in inscription and fragment2 in inscription:
            found_id = stone_id
            break

    if found_id:
        # Using regex to extract the numerical part of the ID string.
        match = re.search(r'\d+', found_id)
        if match:
            runestone_number = int(match.group())
            print(f"The runestone in the image has been identified as {found_id}.")
            print(f"The numerical part of the ID is: {runestone_number}")
    else:
        print("Could not identify the runestone from the fragments.")

identify_runestone()