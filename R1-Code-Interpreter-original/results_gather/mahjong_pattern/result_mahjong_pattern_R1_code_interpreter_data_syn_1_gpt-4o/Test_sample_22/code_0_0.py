from collections import Counter

def determine_result(cards, new_card):
    # Add the new card
    cards.append(new_card)
    
    # Check for "Peng"
    card_count = Counter(cards)
    if card_count[new_card] >= 3:
        return 1  # Peng
    
    # Check for "Chi"
    sorted_cards = sorted(set(cards))
    new_card_index = sorted_cards.index(new_card)
    
    # Check for consecutive sequence
    if new_card_index > 0 and new_card_index < len(sorted_cards) - 1:
        if (ord(sorted_cards[new_card_index - 1]) == ord(new_card) - 1 and
            ord(sorted_cards[new_card_index + 1]) == ord(new_card) + 1):
            return 2  # Chi
    
    return 0  # Pass

# Initial cards
cards = list("EGQTSOAEGBQLU")

# Round 1
cards.remove('T')
determine_result(cards, 'Y')

# Round 2
cards.remove('E')
determine_result(cards, 'E')

# Round 3
cards.remove('Y')
determine_result(cards, 'Q')

# Round 4
cards.remove('Q')
determine_result(cards, 'G')

# Round 5
result = determine_result(cards, 'Q')

print(result)