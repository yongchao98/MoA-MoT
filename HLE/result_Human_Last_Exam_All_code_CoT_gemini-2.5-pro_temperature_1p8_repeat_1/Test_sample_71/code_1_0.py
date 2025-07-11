def solve_old_russian_stress():
    """
    This function determines the stressed syllable in a list of Old Russian phrases
    based on a derived set of linguistic rules.
    """
    phrases = [
        'i ne znali',
        'i povelo že',
        'ne vymyla že',
        'ponesla',
        'vyvela že',
        'i unesli'
    ]
    
    # vowels are a, e, i, o, u, y
    vowels = "aeiouy"

    # Store the results for each phrase
    results = []
    
    # 1. 'i ne znali'
    # The verb root is 'znal', which has fixed root stress (Class R).
    # The stress falls on the root 'zna'.
    # Phrase: i(1) ne(2) zna(3)-li(4). The stressed syllable is the 3rd.
    phrase = phrases[0]
    stress_pos = 3
    results.append(stress_pos)
    print(f"'{phrase}': The verb root 'znal' has fixed stress on the root. The stressed syllable is the {stress_pos}rd one ('zna').")
    
    # 2. 'i povelo že'
    # The verb root is 'vel', which has mobile stress (Class M).
    # With particle 'že' and no preceding 'ne', 'že' takes the stress.
    # Phrase: i(1) po(2)-ve(3)-lo(4) že(5). The stressed syllable is the 5th.
    phrase = phrases[1]
    stress_pos = 5
    results.append(stress_pos)
    print(f"'{phrase}': The verb root 'vel' has mobile stress. The particle 'že' takes the stress. The stressed syllable is the {stress_pos}th one ('že').")

    # 3. 'ne vymyla že'
    # The verb has the prefix 'vy-', which always takes the stress (Class V).
    # Phrase: ne(1) vy(2)-my(3)-la(4) že(5). The stressed syllable is the 2nd.
    phrase = phrases[2]
    stress_pos = 2
    results.append(stress_pos)
    print(f"'{phrase}': The prefix 'vy-' always takes the stress. The stressed syllable is the {stress_pos}nd one ('vy').")
    
    # 4. 'ponesla'
    # The verb root is 'nes', which has mobile stress (Class M).
    # With no relevant particles, the stress falls on the ending.
    # Phrase: po(1)-nes(2)-la(3). The stressed syllable is the 3rd.
    phrase = phrases[3]
    stress_pos = 3
    results.append(stress_pos)
    print(f"'{phrase}': The verb root 'nes' has mobile stress. With no particles, the ending is stressed. The stressed syllable is the {stress_pos}rd one ('la').")
    
    # 5. 'vyvela že'
    # The verb has the prefix 'vy-', which always takes the stress (Class V).
    # Phrase: vy(1)-ve(2)-la(3) že(4). The stressed syllable is the 1st.
    phrase = phrases[4]
    stress_pos = 1
    results.append(stress_pos)
    print(f"'{phrase}': The prefix 'vy-' always takes the stress. The stressed syllable is the {stress_pos}st one ('vy').")
    
    # 6. 'i unesli'
    # The verb root is 'nes', which has mobile stress (Class M).
    # In the 'i + verb' case, the stress falls on the ending.
    # Phrase: i(1) u(2)-nes(3)-li(4). The stressed syllable is the 4th.
    phrase = phrases[5]
    stress_pos = 4
    results.append(stress_pos)
    print(f"'{phrase}': The verb root 'nes' has mobile stress. With only particle 'i', the ending is stressed. The stressed syllable is the {stress_pos}th one ('li').")
    
    # Combine the results into a single string
    final_answer = "".join(map(str, results))
    print("\nThe final sequence of digits is:")
    print(final_answer)

solve_old_russian_stress()