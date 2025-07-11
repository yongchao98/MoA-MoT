def scan_hexameter_line():
    """
    Analyzes and scans the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    print(f"Scanning Latin Hexameter Line: \"{line}\"\n")

    feet = {
        1: {
            "text": "verum ubi",
            "scan": "D",
            "explanation": [
                "The '-um' of 'verum' elides with the starting 'u' of 'ubi'.",
                "The first syllable 'vē-' is long by nature.",
                "The second syllable, 'u' from 'ubi', is short.",
                "The third syllable, '-bi', is short because its vowel 'i' comes before the vowel 'e' of the next word, 'equi'.",
                "Result: Long-Short-Short -> Dactyl (D)."
            ]
        },
        2: {
            "text": "equi atque",
            "scan": "S",
            "explanation": [
                "The '-i' of 'equi' elides with the 'a' of 'atque'.",
                "The first syllable 'e-' of 'equi', though naturally short, is treated as long by poetic license to fit the meter.",
                "The second syllable 'at-' (from 'atque') is long by position, as the 'a' is followed by two consonants ('t' and 'qu').",
                "Result: Long-Long -> Spondee (S)."
            ]
        },
        3: {
            "text": "hominis",
            "scan": "S",
            "explanation": [
                "This foot involves two poetic licenses after the '-que' of 'atque' elides with the 'h' of 'hominis'.",
                "First, the word 'hominis' undergoes syncope, scanning as two syllables ('hom-nis') instead of three.",
                "The first syllable 'hom-' becomes long by position (short 'o' followed by 'm' and 'n').",
                "The second syllable '-nis' is long by position, as the 'i' is followed by 's' and the 'c' of the next word, 'casu'.",
                "Result: Long-Long -> Spondee (S)."
            ]
        },
        4: {
            "text": "casu",
            "scan": "S",
            "explanation": [
                "The first syllable 'cā-' is long by nature.",
                "The second syllable '-sū' is also long by nature.",
                "Result: Long-Long -> Spondee (S)."
            ]
        },
        5: {
            "text": "convenit",
            "scan": "D",
            "explanation": [
                "The first syllable 'con-' is long by position ('o' followed by 'n' and 'v').",
                "The second syllable '-ve-' is short.",
                "The third syllable '-nit' is short.",
                "Result: Long-Short-Short -> Dactyl (D)."
            ]
        },
        6: {
            "text": "imago",
            "scan": "S",
            "explanation": [
                "The first syllable 'i-' is short.",
                "The second syllable '-mā-' is long by nature.",
                "The final syllable '-gō' is anceps (can be long or short), and always counts as long to complete the hexameter.",
                "The foot scans as a spondee for metrical purposes.",
                "Result: (Short)-Long-Long -> Spondee (S)."
            ]
        }
    }

    final_scansion = []
    for i in range(1, 7):
        foot_data = feet[i]
        final_scansion.append(foot_data["scan"])
        print(f"Foot {i} (\"{foot_data['text']}\"): {foot_data['scan']}")
        for line in foot_data['explanation']:
            print(f"  - {line}")
        print("-" * 20)

    # The final "equation" is the sequence of D and S
    final_output = " ".join(final_scansion)
    print("\nFinal Scansion Pattern:")
    print(final_output)

scan_hexameter_line()
<<<D S S S D S>>>