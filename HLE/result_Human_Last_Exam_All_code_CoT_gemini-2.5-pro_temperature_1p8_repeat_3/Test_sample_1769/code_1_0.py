def find_hypothetical_word():
    """
    This function traces the linguistic evolution of the Old English word
    for "sister" to its hypothetical modern form.
    """
    old_english_word = "sweostor"
    norse_equivalent = "systir"
    modern_word = "sister"

    print(f"The Modern English word '{modern_word}' came from the Old Norse '{norse_equivalent}'.")
    print(f"Without Norse influence, the native Old English word '{old_english_word}' would have evolved naturally.")
    print("-" * 30)
    print("Tracing the phonetic evolution of 'sweostor':\n")

    # Step 1: Analyze the 'sweo-' prefix
    print("1. The 'sweo-' prefix evolution:")
    print("   - Compare to Old English 'sweord' which became Modern English 'sword'.")
    print("   - Following this rule: 'sweo-' -> 'swo-'")
    part1 = "swo"

    # Step 2: Analyze the '-stor' suffix
    print("\n2. The '-or' ending evolution:")
    print("   - Compare to Old English 'broþor' which became Modern English 'brother'.")
    print("   - Following this rule, the ending '-or' -> '-er'.")
    print("   - Thus, '-stor' -> '-ster'.")
    part2 = "ster"

    # Step 3: Combine the parts to form the final word
    hypothetical_word = part1 + part2
    print("\n" + "-" * 30)
    print(f"Combining the evolved parts '{part1}' and '{part2}' gives us the final word.")
    print(f"\nThe logical evolution is: {old_english_word} -> {hypothetical_word}")

find_hypothetical_word()