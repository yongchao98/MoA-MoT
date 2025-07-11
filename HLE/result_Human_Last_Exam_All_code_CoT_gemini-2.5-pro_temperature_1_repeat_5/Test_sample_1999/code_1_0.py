import sys

# Plan:
# 1. Determine the order of the five Spice Girls based on the rap bridge in "Wannabe".
# 2. Identify the squares on a chessboard's first rank, starting from White's queenside.
# 3. Map the Spice Girls to these squares in order.
# 4. Find the White Queen's starting square within this sequence.
# 5. Output the nickname prefix of the Spice Girl standing on the Queen's square.

# Step 1: Define the order of the Spice Girls based on the rap references.
# "Em" -> Baby, "G" -> Ginger, "Easy V" -> Posh, "me" -> Sporty, unmentioned -> Scary.
spice_girls_nicknames = ["Baby", "Ginger", "Posh", "Sporty", "Scary"]

# Step 2: Define the relevant chessboard squares and piece names.
# The lineup starts at White's Queenside Rook (a1) and continues along the rank.
files = ['a', 'b', 'c', 'd', 'e']
rank = "1"
squares = [f + rank for f in files]
piece_names = ["Queenside Rook", "Queenside Knight", "Queenside Bishop", "Queen", "King"]

# Step 3 & 4: Map the members to the squares and find the one on the Queen's square.
# The code will iterate through the positions, print the mapping for each,
# and identify the member on the target square.
print("The Spice Girls are placed on the chessboard's first rank as follows:")
target_member = ""
queen_square = "d1"

for i in range(len(squares)):
    square = squares[i]
    member = spice_girls_nicknames[i]
    piece_name = piece_names[i]

    # As requested, we output the mapping for each position.
    print(f"Position {i+1} (Square {square}, White's {piece_name}) = {member} Spice")

    # Check if the current square is the Queen's starting position.
    if square == queen_square:
        target_member = member

# Step 5: Output the final result.
print(f"\nThe White Queen's starting square is {queen_square}.")
print(f"The Spice Girl on this square is {target_member} Spice.")
print(f"The answer is the single word from her nickname that comes before 'Spice'.")

# The final answer is printed in the required format.
sys.stdout.write(f"<<<{target_member}>>>\n")