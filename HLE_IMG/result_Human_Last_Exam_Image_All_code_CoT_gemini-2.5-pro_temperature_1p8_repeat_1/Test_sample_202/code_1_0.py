#
#
# This script identifies the rock specimens based on their visual characteristics.
#
def identify_rocks():
    """
    Identifies and prints the name of each rock specimen shown in the image.
    """
    # A dictionary mapping specimen labels to their identified rock type.
    specimens = {
        'H': 'Scoria',
        'I': 'Breccia',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for label, rock_type in specimens.items():
        print(f"Specimen {label}: {rock_type}")

if __name__ == "__main__":
    identify_rocks()