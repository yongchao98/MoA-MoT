import math

def solve_chair_puzzle():
    """
    Solves the chair puzzle by determining the optimal seating strategy.
    """
    
    # The total number of chairs in the row.
    num_chairs = 20

    # The puzzle's rules lead to a key insight: the number of occupied chairs
    # only increases if a person sits in a chair with no occupied neighbors.
    # To get the maximum number of occupants, we must always follow this strategy.
    
    # This results in a final configuration where no two occupied chairs are adjacent.
    # The most dense pattern of this type is [Occupied, Empty, Occupied, Empty, ...].
    # In this pattern, each person occupies a "block" of 2 chairs (one for them, one empty).
    
    chairs_per_person_in_pattern = 2
    
    # To find the maximum number of people, we divide the total number of chairs
    # by the number of chairs needed for each person in the optimal pattern.
    max_people = num_chairs / chairs_per_person_in_pattern
    
    print("To maximize the number of occupied chairs, people must be seated so that no two are adjacent.")
    print("This corresponds to placing one person on every other chair.")
    print("\nFor a row of 20 chairs, the calculation is:")
    
    # The final equation with each number explicitly shown, as requested.
    print(f"{num_chairs} chairs / {chairs_per_person_in_pattern} chairs per person = {int(max_people)} people")

solve_chair_puzzle()