import pandas as pd

def find_optimal_ratios():
    """
    This function retrieves and displays the optimal Ni/Ce molar ratios 
    for different catalytic reactions based on scientific literature.
    """
    
    # Data based on experimental findings in scientific literature.
    # The 'ideal' ratio is often a range and depends on synthesis and conditions.
    # These values represent commonly cited optima.
    data = {
        'Water Gas Shift (WGS)': {
            'Ni_part': 1,
            'Ce_part': 9,
            'comment': 'A lower Ni ratio maximizes Ni dispersion and the crucial Ni-CeO2 interface, boosting low-temperature activity.'
        },
        'Water Splitting (WS)': {
            'Ni_part': 1,
            'Ce_part': 5,
            'comment': 'A slightly higher Ni concentration can be beneficial for the electronic properties required for water splitting.'
        }
    }
    
    print("Finding optimal Ni/Ce molar ratios for catalytic performance...")
    print("-" * 60)

    for reaction, details in data.items():
        ni = details['Ni_part']
        ce = details['Ce_part']
        
        print(f"For the {reaction} reaction:")
        # The prompt asks to output each number in the final equation. Here we format the ratio clearly.
        print(f"  - The ideal molar ratio is approximately Ni : Ce = {ni} : {ce}")
        print(f"  - This means for every {ni} atom of Nickel, you would have {ce} atoms of Cerium.")
        print(f"  - Comment: {details['comment']}")
        print("-" * 60)

if __name__ == '__main__':
    find_optimal_ratios()