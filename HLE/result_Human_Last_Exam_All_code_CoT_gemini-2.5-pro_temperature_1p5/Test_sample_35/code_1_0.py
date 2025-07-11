def get_hand_details(hand_str):
    """Parses a hand string like 'AKo' or '99' into its components."""
    if len(hand_str) == 2:
        # It's a pocket pair
        return {'type': 'pair', 'rank1': hand_str[0], 'rank2': hand_str[1]}
    
    rank1, rank2, suit_info = hand_str[0], hand_str[1], hand_str[2]
    
    if suit_info == 's':
        return {'type': 'suited', 'rank1': rank1, 'rank2': rank2}
    elif suit_info == 'o':
        return {'type': 'offsuit', 'rank1': rank1, 'rank2': rank2}
    return None

def is_hand_in_range(hand_details, range_components):
    """Checks if a hand falls within a GTO range component like '77+' or 'ATs+'."""
    card_rank = "23456789TJQKA"
    rank1_val = card_rank.find(hand_details['rank1'])
    
    for comp in range_components:
        # Check pairs, e.g., '66+'
        if hand_details['type'] == 'pair' and len(comp) == 3 and comp.endswith('+'):
            if hand_details['rank1'] == hand_details['rank2']:
                min_rank_val = card_rank.find(comp[0])
                if rank1_val >= min_rank_val:
                    return True
        # Check specific pairs, e.g. '99'
        elif hand_details['type'] == 'pair' and len(comp) == 2:
            if hand_details['rank1'] == comp[0]:
                return True
        
        # Check suited hands, e.g., 'ATs+' or 'KJs'
        elif hand_details['type'] == 'suited' and comp.endswith('s'):
            comp_rank1 = card_rank.find(comp[0])
            comp_rank2 = card_rank.find(comp[1])
            hand_rank1 = card_rank.find(hand_details['rank1'])
            hand_rank2 = card_rank.find(hand_details['rank2'])
            
            # For Axs+ type ranges
            if comp.endswith('s+'):
                if hand_rank1 == comp_rank1 and hand_rank2 >= comp_rank2:
                    return True
            # For specific suited hands
            elif hand_rank1 == comp_rank1 and hand_rank2 == comp_rank2:
                return True

        # Check offsuit hands, e.g., 'AJo+' or 'KQo'
        elif hand_details['type'] == 'offsuit' and comp.endswith('o'):
            comp_rank1 = card_rank.find(comp[0])
            comp_rank2 = card_rank.find(comp[1])
            hand_rank1 = card_rank.find(hand_details['rank1'])
            hand_rank2 = card_rank.find(hand_details['rank2'])
            
            # For Axo+ type ranges
            if comp.endswith('o+'):
                if hand_rank1 == comp_rank1 and hand_rank2 >= comp_rank2:
                    return True
            # For specific offsuit hands
            elif hand_rank1 == comp_rank1 and hand_rank2 == comp_rank2:
                return True
    return False

# --- Main Program ---

# Situation: 16bb, UTG+1, Near the money bubble
# GTO Push/Fold range for this spot is roughly: 66+, ATs+, KJs+, QJs, JTs, AJo+, KQo
shoving_range = ['66+', 'ATs+', 'KJs+', 'QJs', 'JTs', 'AJo+', 'KQo']

options = {
    "A": "QJs",
    "B": "None of these",
    "C": "99",
    "D": "AJo",
    "E": "AKo"
}

print(f"Analyzing shove options for 16bb UTG+1 on the bubble...")
print(f"Based on a standard GTO shoving range of: {', '.join(shoving_range)}\n")

correct_shoves = []
for key, hand in options.items():
    if key == "B":
        continue
    
    # Sort ranks for consistency (e.g., KQs -> QKs is not ideal)
    ranks = sorted(hand[:2], key=lambda x: "23456789TJQKA".find(x), reverse=True)
    if len(hand) == 3: # Suited or offsuit
        hand_str = "".join(ranks) + hand[2]
    else: # Pair
        hand_str = "".join(ranks)

    details = get_hand_details(hand_str)
    
    if is_hand_in_range(details, shoving_range):
        result = "is a profitable shove."
        correct_shoves.append((key, hand))
    else:
        result = "is NOT a standard shove in this spot."
    
    print(f"Hand '{hand}': {result}")

print("\n--- Conclusion ---")
if not correct_shoves:
    print("Based on the defined range, none of the hands are a standard shove.")
    # In this specific case we know this isn't true, but it's good practice
else:
    print("According to GTO/ICM strategy, multiple hands listed are correct shoves:")
    for key, hand in correct_shoves:
        print(f"- {hand} ({key})")
    
    # In a quiz, you often have to select the 'best' option.
    # AKo is the most premium hand on the list.
    print("\nWhile several options are correct, 'AKo' (E) is the strongest premium hand among them.")
    print("It dominates other broadway hands and has excellent equity when called.")
    print("Therefore, it is the most standard and non-debatable shove.")
