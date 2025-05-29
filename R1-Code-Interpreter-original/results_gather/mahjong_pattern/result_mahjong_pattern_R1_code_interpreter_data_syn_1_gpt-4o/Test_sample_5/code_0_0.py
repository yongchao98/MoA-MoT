# Initial cards
cards = list("YQILCUCZEAXNK")

# Round 1
cards.append("H")  # Add H
cards.remove("E")  # Remove E

# Check for "Peng" in Round 1
def check_peng(cards, new_card):
    card_count = {}
    for card in cards:
        if card in card_count:
            card_count[card] += 1
        else:
            card_count[card] = 1
    # Check if there are two identical cards and the new card is the same
    for card, count in card_count.items():
        if count >= 2 and card == new_card:
            return 1  # Peng
    return 0  # Pass

# Round 2
cards.append("C")  # Add C

# Determine the result for the final round
result = check_peng(cards, "C")
print(result)