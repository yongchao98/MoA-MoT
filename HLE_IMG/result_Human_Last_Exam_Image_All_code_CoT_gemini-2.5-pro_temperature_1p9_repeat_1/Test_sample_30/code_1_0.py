def check_solvability(config_name, config_data):
    """
    Checks if a configuration is solvable based on a parity invariant.
    
    A necessary condition for solvability is that the number of black knights
    on black-colored squares (Nb_b) and the number of white knights on
    black-colored squares (Nw_b) must have the same parity.
    
    Args:
        config_name (str): The label of the configuration (e.g., 'A').
        config_data (dict): A dictionary with the locations of black and white knights.

    Returns:
        bool: True if the configuration might be solvable, False otherwise.
    """
    
    # A square (r, c) is 'black-colored' if (r+c) is even.
    
    num_black_on_black_sq = 0
    for r, c in config_data['black']:
        if (r + c) % 2 == 0:
            num_black_on_black_sq += 1
            
    num_white_on_black_sq = 0
    for r, c in config_data['white']:
        if (r + c) % 2 == 0:
            num_white_on_black_sq += 1
            
    is_solvable = (num_black_on_black_sq % 2 == num_white_on_black_sq % 2)
    
    print(f"Configuration {config_name}:")
    print(f"  Number of black knights on black-colored squares: {num_black_on_black_sq}")
    print(f"  Number of white knights on black-colored squares: {num_white_on_black_sq}")
    print(f"  Parity of black knights on black squares: {num_black_on_black_sq % 2}")
    print(f"  Parity of white knights on black squares: {num_white_on_black_sq % 2}")
    
    if is_solvable:
        print(f"  Parities are the same. Configuration {config_name} is potentially SOLVABLE.")
    else:
        print(f"  Parities are different. Configuration {config_name} is UNSOLVABLE.")
    print("-" * 20)
    
    return is_solvable


def main():
    """
    Defines the initial configurations from the image and checks their solvability.
    Coordinates are 0-indexed from the top-left.
    """
    configs = {
        'A': {
            'black': [(0, 0), (1, 0), (2, 0), (3, 0)],
            'white': [(0, 2), (1, 2), (2, 2), (3, 2)]
        },
        'B': {
            'black': [(0, 0), (2, 0), (2, 2)],
            'white': [(1, 1), (3, 1), (3, 2)]
        },
        'C': {
            'black': [(0, 2), (1, 2)],
            'white': [(0, 0), (2, 1)]
        },
        'D': {
            'black': [(1, 1), (3, 1)],
            'white': [(0, 0), (2, 1)]
        },
        'E': {
            'black': [(0, 0), (1, 0), (1, 1)],
            'white': [(0, 1), (0, 2), (1, 2)]
        }
    }
    
    solvable_configs = []
    for name, data in configs.items():
        if check_solvability(name, data):
            solvable_configs.append(name)
            
    print(f"\nBased on the parity invariant, the solvable configurations are: {', '.join(solvable_configs)}")


if __name__ == "__main__":
    main()
