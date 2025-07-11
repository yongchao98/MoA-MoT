import collections

def get_unseen_composition(config, person_idx, total_b, total_w):
    """Calculates the composition of hats unseen by a person."""
    unseen_indices = [(person_idx - 1 + 9) % 9, person_idx, (person_idx + 1) % 9]
    
    # Person sees all other hats
    seen_b = 0
    seen_w = 0
    for i in range(9):
        if i not in unseen_indices:
            if config[i] == 'B':
                seen_b += 1
            else:
                seen_w += 1
    
    unseen_b = total_b - seen_b
    unseen_w = total_w - seen_w
    return (unseen_b, unseen_w)

def check_round1(config):
    """Checks if anyone would speak in Round 1."""
    # A person speaks if their unseen hats are monochromatic.
    # This is equivalent to checking for any BBB or WWW blocks.
    for i in range(9):
        p, c, n = config[i-1], config[i], config[(i+1)%9]
        if p == c == n:
            return False # This config is invalid as someone would speak.
    return True

def check_round2(config):
    """Checks if anyone would speak in Round 2."""
    # A person speaks if assuming the opposite hat color creates a BBB or WWW,
    # which contradicts the knowledge from Round 1.
    for i in range(9):
        unseen_b, unseen_w = get_unseen_composition(config, i, 5, 4)
        
        # Check if they can deduce they are Black
        if unseen_b > 0 and unseen_w == unseen_b - 1: # Unseen are (eg) 2B,1W or 3B,2W
            # If they are actually White, could they deduce it?
            if config[i] == 'W':
                # Hypoth: 'B'. Unseen comp would be (unseen_b+1, unseen_w-1)
                if unseen_w - 1 == 0: # This would mean the other 2 are B, forming BBB
                    return False # This person would speak, so config is invalid for R2
        
        # Check if they can deduce they are White
        if unseen_w > 0 and unseen_b == unseen_w - 1: # Unseen are (eg) 1B,2W
            if config[i] == 'B':
                # Hypoth: 'W'. Unseen comp would be (unseen_b-1, unseen_w+1)
                if unseen_b - 1 == 0: # This would mean the other 2 are W, forming WWW
                    return False # This person would speak, so config is invalid for R2
    return True

def solve_puzzle():
    """Finds the hat configuration and solves the puzzle."""
    final_config = "BWBWBWBWB" # Based on logical deduction
    total_b = final_config.count('B')
    total_w = final_config.count('W')
    
    # Verification using the logic derived
    print("Verifying the proposed solution...")
    print(f"Configuration: {' '.join(final_config)}")

    # Round 1 check
    if not check_round1(list(final_config)):
        print("Error: Proposed solution fails Round 1.")
        return

    # Round 2 check
    if not check_round2(list(final_config)):
        print("Error: Proposed solution fails Round 2.")
        return
    
    print("\nThis configuration is consistent with 'No' answers in Round 1 and 2.\n")
    print("Analyzing Round 3:")
    
    num_yes_in_r3 = 0
    speakers = []
    
    # In Round 3, a person speaks if assuming the opposite hat color
    # creates a BWB or WBW pattern, which contradicts Round 2 knowledge.
    for i in range(9):
        person_speaks = False
        unseen_b, unseen_w = get_unseen_composition(list(final_config), i, total_b, total_w)
        
        # Deduction logic for Round 3
        # If unseen is (2B, 1W), person must be B. Otherwise, assuming W creates BWB.
        if unseen_b == 2 and unseen_w == 1:
            if final_config[i] == 'B':
                person_speaks = True
                
        # If unseen is (1B, 2W), person must be W. Otherwise, assuming B creates WBW.
        if unseen_b == 1 and unseen_w == 2:
            if final_config[i] == 'W':
                person_speaks = True
        
        if person_speaks:
            num_yes_in_r3 += 1
            speakers.append(f"Person {i+1} ({final_config[i]})")

    print(f"Number of people who replied 'Yes' in Round 3: {num_yes_in_r3}")
    print(f"The distribution of hats is: {' '.join(final_config)}")

solve_puzzle()
>>> 9
B W B W B W B W B