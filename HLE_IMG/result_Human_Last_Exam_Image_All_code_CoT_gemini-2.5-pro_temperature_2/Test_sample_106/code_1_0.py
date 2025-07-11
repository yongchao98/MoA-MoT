def solve_rearrangement():
    """
    This function determines and prints the identity of the substituents at
    positions 1, 2, 3, 4, and 5 in the final product of the rearrangement reaction.
    """
    substituents = {
        '1': 'CH3',
        '2': 'H',
        '3': 'CH3',
        '4': 'CH3',
        '5': 'H'
    }
    
    output = []
    for i in range(1, 6):
        key = str(i)
        value = substituents[key]
        output.append(f"{key} = {value}")
        
    print(", ".join(output))

solve_rearrangement()