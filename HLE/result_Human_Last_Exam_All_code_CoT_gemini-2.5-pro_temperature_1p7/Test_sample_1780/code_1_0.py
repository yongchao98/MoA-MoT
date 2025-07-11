import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters "L", "N", "S", and "W" can be arranged
    based on a specific tip-connection rule.
    """
    
    # Step 1: Define the properties of each letter's tips.
    # The tuple represents (Leftmost Tip, Rightmost Tip).
    letter_tips = {
        'L': ('bottom', 'bottom'),
        'N': ('top', 'bottom'),
        'S': ('top', 'bottom'),
        'W': ('top', 'top')
    }
    
    letters = ['L', 'N', 'S', 'W']
    
    # Step 2: Generate all possible arrangements (permutations).
    all_arrangements = itertools.permutations(letters)
    
    valid_arrangement_count = 0
    
    # Step 3: Check each arrangement against the connection rule.
    for arr in all_arrangements:
        is_valid = True
        # Check connections between adjacent letters.
        for i in range(len(arr) - 1):
            letter1 = arr[i]
            letter2 = arr[i+1]
            
            right_tip1 = letter_tips[letter1][1]
            left_tip2 = letter_tips[letter2][0]
            
            # Step 4: The rule, derived from the examples, is that tips must be different.
            if right_tip1 == left_tip2:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangement_count += 1
            
    # Step 5: Output the result as an equation showing the sum.
    if valid_arrangement_count > 0:
        sum_parts = " + ".join(['1'] * valid_arrangement_count)
        print(f"{sum_parts} = {valid_arrangement_count}")
    else:
        print("0")

solve_letter_arrangement()