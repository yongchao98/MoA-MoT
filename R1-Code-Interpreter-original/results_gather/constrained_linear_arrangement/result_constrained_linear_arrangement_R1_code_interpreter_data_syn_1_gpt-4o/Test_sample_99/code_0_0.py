# Deduce opponent's cards based on the given rounds and results
def deduce_opponent_cards():
    # Initialize the opponent's cards
    opponent_cards = [''] * 5
    
    # Round 1: H, loss
    opponent_cards[0] = 'G'  # Opponent played G to defeat H
    
    # Round 2: F, draw
    opponent_cards[1] = 'F'  # Opponent played F for a draw
    
    # Round 3: C, win
    opponent_cards[2] = 'D'  # Opponent played D and lost to C
    
    # Round 4: A, win
    opponent_cards[3] = 'B'  # Opponent played B and lost to A
    
    # Round 5: E, win
    opponent_cards[4] = 'H'  # Opponent played H and lost to E
    
    return opponent_cards

# Print the deduced opponent's cards
print(deduce_opponent_cards())