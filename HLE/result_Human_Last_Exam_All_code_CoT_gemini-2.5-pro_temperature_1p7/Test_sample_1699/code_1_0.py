import collections

def solve_go_problem():
    """
    Analyzes the Go board position and determines the best move for Black to capture all White stones.
    """
    # Initial configuration of stones
    # Each piece is (row, column), where row is 1-19 (top-bottom) and column is 1-19 (right-left).
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
    
    print("Analyzing the Go board to find the move for Black to capture all White stones.")
    print("-" * 30)

    # --- Step 1: Analyze the vital point ---
    # The white stones are in a loose arrangement, trying to form a living group.
    # To live, a group needs to make two separate "eyes" (protected empty points).
    # Black's goal is to prevent this by taking the most critical point in White's potential eye space.
    # The point (2, 4) is a nexus, connecting multiple white stones and groups.
    
    best_move = (2, 4)
    print(f"Step 1: Black plays at the vital point {best_move}.")
    print("This move is strategically powerful because it splits the White stones and immediately creates a threat.")
    
    # Update board state
    black_stones.add(best_move)
    
    # --- Step 2: Analyze the immediate consequence ---
    # Check the liberties of the white stone at (2, 5) after Black's move.
    # A liberty is an adjacent empty point.
    
    w_stone_to_check = (2, 5)
    neighbors = {
        (w_stone_to_check[0] - 1, w_stone_to_check[1]),  # up
        (w_stone_to_check[0] + 1, w_stone_to_check[1]),  # down
        (w_stone_to_check[0], w_stone_to_check[1] - 1),  # right
        (w_stone_to_check[0], w_stone_to_check[1] + 1)   # left
    }
    
    liberties = []
    for neighbor in neighbors:
        if neighbor not in black_stones and neighbor not in white_stones:
            liberties.append(neighbor)
            
    print(f"\nAfter Black's move, let's check the white stone at {w_stone_to_check}.")
    print(f"Its neighbors are B(3,5), B(2,6), B(2,4) and the empty point (1,5).")
    print(f"This puts the white stone at {w_stone_to_check} in 'atari'. It has only one liberty left at {liberties[0]}.")

    # --- Step 3: White's forced move ---
    white_response = (1, 5)
    print(f"\nStep 2: White must play at {white_response} to save the stone at {w_stone_to_check}.")
    print("If White plays anywhere else, Black will capture the stone on the next turn.")
    white_stones.add(white_response)
    
    # --- Step 4: Black's follow-up ---
    black_follow_up = (1, 3)
    print(f"\nStep 3: Black maintains the advantage and plays at {black_follow_up}.")
    print("This move further reduces the space White has to form eyes and puts another white group under pressure.")
    black_stones.add(black_follow_up)
    
    # --- Step 5: Conclusion ---
    print("\nConclusion:")
    print("After this sequence, the White group is split and confined.")
    print("The main white group now has a shape known as a 'bulky five', but Black can occupy its vital point.")
    print("White does not have enough space ('liberties' or 'me' in Japanese) to create the two separate eyes required for life.")
    print("Therefore, White's entire group will be captured.")
    
    print("-" * 30)
    print(f"The first move that makes this possible is Black playing at (2, 4).")
    
if __name__ == '__main__':
    solve_go_problem()