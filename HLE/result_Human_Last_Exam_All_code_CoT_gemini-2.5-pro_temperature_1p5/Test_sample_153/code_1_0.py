import re

def solve_riddle():
    """
    Solves the Cold War riddle by interpreting "Кома" as a numerical code.
    """
    # Step 1: Define the 32-letter Russian alphabet (omitting Ё) and create a mapping.
    alphabet = "АБВГДЕЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ"
    letter_to_pos = {letter: i + 1 for i, letter in enumerate(alphabet)}

    # Step 2: Define the clue word and the location choices.
    clue_word = "КОМА"
    locations = {
        "A": "Калининградская область",
        "B": "Пермский край",
        "C": "Таймырский Долгано-Ненецкий район",
        "D": "Чукотский автономный округ",
        "E": "Республика Адыгея"
    }

    # Step 3: Calculate the target sum from the clue word.
    clue_letters = list(clue_word)
    clue_values = [letter_to_pos[l] for l in clue_letters]
    target_sum = sum(clue_values)

    print(f"Analyzing the clue word: '{clue_word}'")
    clue_calc_str = ' + '.join(map(str, clue_values))
    print(f"Numerical value: {clue_calc_str} = {target_sum}\n")
    print("Searching for a location with a matching value...")
    print("-" * 40)

    # Step 4 & 5: Iterate through locations, get initials, and check their sum.
    for key, name in locations.items():
        # Get the first letter of each word.
        initials = re.findall(r'\b\w', name.upper())
        
        try:
            initial_values = [letter_to_pos[i] for i in initials]
            current_sum = sum(initial_values)
            
            print(f"Location {key}: {name}")
            print(f"Initials: {''.join(initials)}")
            calc_str = ' + '.join(map(str, initial_values))
            print(f"Sum: {calc_str} = {current_sum}")

            if current_sum == target_sum:
                print("\n" + "="*40)
                print("!!! MATCH FOUND !!!")
                print(f"The destination is: {name}")
                print("The decoded instruction establishes the following equality:")
                print(f"Equation: {' + '.join(initials)} = {' + '.join(clue_letters)}")
                print(f"Numerical representation:")
                final_equation = f"{calc_str} = {clue_calc_str}"
                print(final_equation)
                final_result = f"{current_sum} = {target_sum}"
                print(final_result)
                print("="*40)
                return key

        except KeyError as e:
            print(f"Location {key}: {name} -> Could not process initial {e}, skipping.")
        finally:
            print("-" * 40)
            
solve_riddle()