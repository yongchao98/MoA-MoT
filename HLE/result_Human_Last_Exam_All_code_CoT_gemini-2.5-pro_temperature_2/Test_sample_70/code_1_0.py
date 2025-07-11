import sys
from io import StringIO

def get_rank(flag_matrix):
    """
    Calculates the maximal rank of a flag matrix.
    Under the "maximal rank" assumption, the rank is the number of
    unique row patterns in the flag's design.
    """
    # Convert each row to a tuple to make it hashable
    # Then use a set to find the unique rows.
    # The size of the set is the number of unique patterns, which is the rank.
    unique_rows = {tuple(row) for row in flag_matrix}
    return len(unique_rows)

def solve():
    """
    Finds African flags with the same matrix rank as Denmark's flag.
    """
    # Define flag structures as matrices of characters representing colors.
    # The dimensions are chosen to capture the essential color patterns.
    flags = {
        # --- Target Flag ---
        "Denmark": [
            ['R', 'W', 'R'],
            ['W', 'W', 'W'],
            ['R', 'W', 'R']
        ],

        # --- African Flags ---
        "Benin": [  # Vertical green, horizontal yellow/red
            ['G', 'Y', 'Y'],
            ['G', 'R', 'R']
        ],
        "Madagascar": [ # Vertical white, horizontal red/green
            ['W', 'R', 'R'],
            ['W', 'G', 'G']
        ],
        "Nigeria": [ # Vertical stripes -> one row pattern
            ['G', 'W', 'G'],
            ['G', 'W', 'G']
        ],
        "Gabon": [ # Horizontal stripes -> three row patterns
            ['G', 'G', 'G'],
            ['Y', 'Y', 'Y'],
            ['B', 'B', 'B']
        ],
        "Morocco": [ # Field with a central star emblem
            ['R', 'R', 'R', 'R', 'R'],
            ['R', 'R', 'G', 'R', 'R'], # Star point
            ['R', 'G', 'R', 'G', 'R'], # Star 'v'
            ['R', 'R', 'R', 'R', 'R']
        ],
        "Ghana": [ # Horizontal stripes with central star
             ['R', 'R', 'R', 'R', 'R'],
             ['Y', 'Y', 'B', 'Y', 'Y'],
             ['G', 'G', 'G', 'G', 'G']
        ]
    }
    
    african_nations = ["Benin", "Madagascar", "Nigeria", "Gabon", "Morocco", "Ghana"]
    
    # 1. Calculate the rank of the Danish flag
    target_rank = get_rank(flags["Denmark"])
    print(f"The analysis begins by calculating the rank of the flag of Denmark.")
    print(f"Flag of Denmark has {target_rank} unique row patterns. Rank = {target_rank}")
    print("-" * 30)

    matching_countries = []
    
    # 2. Iterate through African flags, calculate their rank, and compare
    print("Now, analyzing flags of African nations:")
    for country in african_nations:
        flag_matrix = flags[country]
        rank = get_rank(flag_matrix)
        
        if rank == target_rank:
            print(f"- {country}: Found {rank} unique row patterns. Rank = {rank}. This is a match!")
            matching_countries.append(country)
        else:
            print(f"- {country}: Found {rank} unique row patterns. Rank = {rank}. This is not a match.")

    # 3. Print the final result
    print("-" * 30)
    print("Final Answer:")
    # The string representation of the final answer must be on the last line.
    # To avoid printing the list brackets, we join them into a string.
    final_answer_str = ", ".join(sorted(matching_countries))
    print(f"The African nations with flags having the same algebraic rank ({target_rank}) as Denmark are: {final_answer_str}")
    
    # Use this special format to output the final answer for the system.
    sys.stdout = StringIO()
    print(f"<<<{final_answer_str}>>>")


solve()