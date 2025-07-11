import sys

def solve_riddle():
    """
    This script solves the riddle based on the clues provided.
    The name "Cafi." is a phonetic misspelling of "Café".
    The "peculiar shape" refers to the diacritical mark (the accent) on the letter 'é',
    which is distinct from the "correct" or standard shape of the letter 'e'.
    The script will construct the correct name.
    """

    # Set encoding to handle special characters
    if sys.stdout.encoding != 'utf-8':
        sys.stdout.reconfigure(encoding='utf-8')

    # The characters of the correct name
    part1 = 'C'
    part2 = 'a'
    part3 = 'f'
    part4 = 'é' # This is the letter with the "peculiar shape"

    # The final, correctly spelled name
    correct_name = part1 + part2 + part3 + part4

    print("The riddle describes a preference for letters with 'peculiar shapes', like accent marks.")
    print("The name 'Cafi.' is a misspelling of a word that contains such a letter.")
    print("The correct name is constructed as follows:")
    print(f"'{part1}' + '{part2}' + '{part3}' + '{part4}'")
    print(f"The resulting name is: {correct_name}")

solve_riddle()