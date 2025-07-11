import collections

def solve_flag_rank_problem():
    """
    This script identifies African national flags that share the same linear algebraic rank as the Danish flag.
    The rank is determined by the flag's geometric pattern, assuming maximal rank where different colors have different values.
    """

    # Step 1: Define the rank of the Danish flag.
    # The Danish flag's cross design creates two unique row patterns (a field row and a crossbar row).
    # Thus, its rank is 2.
    denmark_rank = 2
    print("The flag of Denmark features a Nordic cross, a design that results in 2 unique row patterns.")
    print(f"Therefore, the target rank is {denmark_rank}.\n")

    # Step 2: Define the structure of various African flags.
    # The 'pattern' key describes the basic layout.
    # 'colors' represents the distinct color fields or stripes.
    # 'emblem' indicates the presence of a symbol that complicates the base pattern.
    african_flags = {
        'Nigeria': {'pattern': 'vertical_stripes', 'colors': ['green', 'white', 'green'], 'emblem': False},
        'Benin': {'pattern': 'hoist_fly', 'fly_stripes': 2, 'emblem': False},
        'Madagascar': {'pattern': 'hoist_fly', 'fly_stripes': 2, 'emblem': False},
        'Botswana': {'pattern': 'horizontal_stripes_fimbriated', 'unique_color_rows': 3, 'emblem': False},
        'Cameroon': {'pattern': 'vertical_stripes', 'colors': ['green', 'red', 'yellow'], 'emblem': True},
        'Ghana': {'pattern': 'horizontal_stripes', 'colors': ['red', 'yellow', 'green'], 'emblem': True},
        'Mali': {'pattern': 'vertical_stripes', 'colors': ['green', 'yellow', 'red'], 'emblem': False},
        'Angola': {'pattern': 'horizontal_stripes', 'colors': ['red', 'black'], 'emblem': True},
        'Tanzania': {'pattern': 'diagonal', 'colors': 4, 'emblem': False},
    }

    def get_rank(flag_name, data):
        """Calculates the rank based on the flag's structural properties."""
        pattern = data['pattern']
        has_emblem = data['emblem']
        rank = 0
        explanation = ""

        if pattern == 'vertical_stripes':
            # Rank is the number of unique vertical stripe colors.
            # An emblem on one stripe makes it unique but doesn't change the rank of a tricolor.
            # However, for simplicity here, we consider any emblem as potential added complexity.
            unique_colors = len(set(data['colors']))
            rank = unique_colors
            if has_emblem:
                rank = max(rank, 3) # Emblems on tricolors (Cameroon) result in rank 3.
            explanation = f"Based on {unique_colors} unique vertical stripe colors, its rank is {rank}."

        elif pattern == 'horizontal_stripes':
            unique_colors = len(set(data['colors']))
            rank = unique_colors
            # An emblem adds complexity, creating new row patterns.
            if has_emblem:
                rank += 1
            explanation = f"Based on {unique_colors} horizontal stripe colors{' with an emblem' if has_emblem else ''}, its rank is {rank}."

        elif pattern == 'horizontal_stripes_fimbriated':
            # Fimbriation (bordering) adds unique rows.
            rank = data['unique_color_rows']
            explanation = f"Based on a field, stripe, and fimbriation, it has {rank} unique row patterns."

        elif pattern == 'hoist_fly':
            # The number of horizontal stripes on the fly determines the number of unique row patterns.
            rank = data['fly_stripes']
            explanation = f"Based on its hoist-fly design with {rank} horizontal stripes, its rank is {rank}."

        elif pattern == 'diagonal':
            # Diagonal designs create many unique rows, resulting in a high rank.
            rank = 5  # Placeholder for "high rank"
            explanation = f"Based on its complex diagonal design, its rank is high (e.g., >{rank-1})."
        
        print(f"Analyzing {flag_name}: {explanation}")
        return rank

    print("--- Analyzing African Flags ---")
    matching_flags = []
    for country, data in african_flags.items():
        rank = get_rank(country, data)
        if rank == denmark_rank:
            matching_flags.append(country)
            print(f"-> Match found: Rank is {rank}.\n")
        else:
            print(f"-> No match: Rank is {rank}.\n")

    # Step 3: Print the final result.
    print("--- Result ---")
    if matching_flags:
        print("The following African nations have flags with the same linear algebraic rank as Denmark (Rank 2):")
        for country in sorted(matching_flags):
            print(f"- {country}")
    else:
        print("No African flags with a rank of 2 were found in the analyzed list.")
    
    # Final answer in the specified format
    final_answer = ", ".join(sorted(matching_flags))
    print(f"<<<{final_answer}>>>")

solve_flag_rank_problem()