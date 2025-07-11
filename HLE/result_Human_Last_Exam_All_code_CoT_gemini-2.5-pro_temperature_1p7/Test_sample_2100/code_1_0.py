def solve_riddle():
    # The man in question is Nikolai Rezanov, the subject of Andrei Voznesensky's rock opera "Juno and Avos".
    # All the clues from the prompt point to him.

    last_name = "Rezanov"
    year_of_death = 1807
    year_of_birth = 1764

    print(f"The poet Andrei Voznesensky wrote about this man in the famous rock opera 'Juno and Avos'.")
    print(f"The man's grave was lost, and his story was revived in the public consciousness during the late 1980s.")
    print(f"The man's last name is: {last_name}")

    # To fulfill the request for an equation, we can calculate his age at death.
    age_at_death = year_of_death - year_of_birth
    
    print("\nHere is an equation related to his life, showing his age at death:")
    print(f"{year_of_death} - {year_of_birth} = {age_at_death}")

solve_riddle()