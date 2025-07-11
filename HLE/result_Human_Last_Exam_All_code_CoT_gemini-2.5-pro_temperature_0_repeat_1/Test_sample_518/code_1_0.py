def solve_riddle():
    """
    Solves the Latin riddle by constructing a word from the first syllables
    of the names mentioned.
    """
    # The four names mentioned in the riddle
    names = ["Penelope", "Dido", "Cadmus", "Remus"]
    
    # Extract the first syllable (first two letters) from each name
    s1 = names[0][:2]
    s2 = names[1][:2]
    s3 = names[2][:2]
    s4 = names[3][:2]
    
    syllables = [s1, s2, s3, s4]
    
    # The final word is the concatenation of the syllables
    punishment = "".join(syllables)
    
    print("The riddle asks us to form a word from the first syllables of four names:")
    print(f"1. The first from 'Penelopes' is '{s1}'")
    print(f"2. The first from 'Didonis' is '{s2}'")
    print(f"3. The first from 'Cadmi' is '{s3}'")
    print(f"4. The first from 'Remi' is '{s4}'")
    print("\nCombining these syllables creates the word for the punishment:")
    
    # Show the equation for forming the word
    print(f"'{s1}' + '{s2}' + '{s3}' + '{s4}' = '{punishment.upper()}'")
    
    print(f"\nThe punishment for stealing from the garden is to be subjected to '{punishment}'.")
    print("This is a Latin verb with a crude meaning, implying a harsh physical punishment.")

solve_riddle()