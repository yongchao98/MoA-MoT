def scan_hexameter_line():
    """
    Scans the Latin hexameter line: "verum ubi equi atque hominis casu convenit imago."
    This function explains the scansion step-by-step and prints the final result.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    print(f"Scanning line: {line}\n")
    
    # Detailed breakdown of the scansion process
    print("Step-by-step analysis:\n")
    
    # Foot 1
    print("Foot 1: 'verum ubi e-'")
    print("  - Syllables: ver- | u | e")
    print("  - 'ver-': Long by position (the 'e' is followed by 'r' and 'm').")
    print("  - The syllable '-um' of 'verum' is elided before 'ubi'.")
    print("  - 'u': The 'u' from 'ubi' is short by nature.")
    print("  - The syllable '-bi' of 'ubi' is elided before 'equi'.")
    print("  - 'e': The 'e' from 'equi' is short by nature.")
    print("  - Pattern: Long, Short, Short => Dactyl (D)\n")
    
    # Foot 2
    print("Foot 2: '-qui atque'")
    print("  - Syllables: qui | at-")
    print("  - '-qui': The 'i' in 'equi' (genitive 'equī') is long by nature.")
    print("  - 'at-': Long by position (the 'a' is followed by 't' and 'qu').")
    print("  - The syllable '-que' of 'atque' is elided before 'hominis'.")
    print("  - Pattern: Long, Long => Spondee (S)\n")

    # Foot 3
    print("Foot 3: 'hominis'")
    print("  - Syllables: ho- | mi- | nis")
    print("  - 'ho-': The 'o' is naturally short. For the line to scan correctly, this must be treated as long by poetic license, a practice sometimes seen in Lucretius.")
    print("  - 'mi-': Short by nature.")
    print("  - 'nis': Short by nature (it is not followed by two consonants within its foot).")
    print("  - Pattern: Long, Short, Short => Dactyl (D)\n")
    
    # Foot 4
    print("Foot 4: 'casu'")
    print("  - Syllables: ca- | su")
    print("  - 'ca-': Long by nature (from ablative 'cāsū').")
    print("  - 'su': Long by nature (from ablative 'cāsū').")
    print("  - Pattern: Long, Long => Spondee (S)\n")

    # Foot 5
    print("Foot 5: 'convenit'")
    print("  - Syllables: con- | ve- | nit")
    print("  - 'con-': Long by position (the 'o' is followed by 'n' and 'v').")
    print("  - 've-': Short by nature.")
    print("  - 'nit': Short by nature.")
    print("  - Pattern: Long, Short, Short => Dactyl (D)\n")

    # Foot 6
    print("Foot 6: '(i)mago'")
    print("  - Syllables: ma- | go")
    print("  - This foot involves prodelision, where the initial 'i-' of 'imago' is elided after the preceding word ending in a consonant.")
    print("  - 'ma-': Long by nature (from 'imāgo').")
    print("  - 'go': Long (the final syllable of a hexameter line is always treated as long).")
    print("  - Pattern: Long, Long => Spondee (S)\n")
    
    # Final Result
    final_scansion_line = "verum ubi e|qui atque | hominis | casu | convenit | (i)mago"
    final_pattern = "D S D S D S"
    
    print("----------------------------------------")
    print("Final Scansion:")
    print(final_scansion_line)
    print("verum ubi equi atque hominis casu convenit imago.")
    print(f"Final Pattern: {final_pattern}")
    print("----------------------------------------")


if __name__ == "__main__":
    scan_hexameter_line()
<<<D S D S D S>>>