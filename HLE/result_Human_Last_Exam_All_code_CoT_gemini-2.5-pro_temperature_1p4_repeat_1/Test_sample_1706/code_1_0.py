def get_neighbors(r, c):
    """Gets valid neighbors for a coordinate on a 19x19 board."""
    neighbors = []
    if r > 1: neighbors.append((r - 1, c))
    if r < 19: neighbors.append((r + 1, c))
    if c > 1: neighbors.append((r, c - 1))
    if c < 19: neighbors.append((r, c + 1))
    return neighbors

def find_group(r, c, stones):
    """Finds all stones in a connected group using BFS."""
    if (r, c) not in stones:
        return set(), set()
    
    q = [(r, c)]
    group = set(q)
    visited = set(q)
    
    head = 0
    while head < len(q):
        curr_r, curr_c = q[head]
        head += 1
        
        for neighbor in get_neighbors(curr_r, curr_c):
            if neighbor in stones and neighbor not in visited:
                visited.add(neighbor)
                group.add(neighbor)
                q.append(neighbor)
    return group

def get_liberties(group, black_stones, white_stones):
    """Calculates the liberties for a given group of stones."""
    liberties = set()
    all_stones = black_stones.union(white_stones)
    for r, c in group:
        for neighbor in get_neighbors(r, c):
            if neighbor not in all_stones:
                liberties.add(neighbor)
    return liberties

def analyze_and_print(turn, move, black_stones, white_stones):
    """Analyzes the board state and prints a summary."""
    player = "Black" if turn % 2 != 0 else "White"
    print(f"--- Turn {turn}: {player} plays at {move} ---")
    
    # Find all distinct white groups and their liberties
    w_groups_found = set()
    white_groups_analysis = []
    for r, c in sorted(list(white_stones)):
        stone_tuple = tuple(sorted(list(find_group(r, c, white_stones))))
        if stone_tuple not in w_groups_found:
            w_groups_found.add(stone_tuple)
            group = set(stone_tuple)
            liberties = get_liberties(group, black_stones, white_stones)
            analysis = {
                "group": group,
                "liberties_count": len(liberties),
                "liberties": liberties
            }
            white_groups_analysis.append(analysis)

    print("Current White groups and their liberties:")
    for item in white_groups_analysis:
        group_str = ", ".join(map(str, sorted(list(item['group']))))
        status = "(ATARI)" if item['liberties_count'] == 1 else ""
        print(f"  - Group {{{group_str}}}: {item['liberties_count']} liberties {status}")
    print("")

def solve_go_problem():
    """Simulates the sequence of moves to solve the problem."""
    # Initial setup
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    print("Initial state analysis:")
    analyze_and_print(0, "Start", black_stones, white_stones)

    # --- Move Sequence ---
    # 1. B:(2,4)
    move1_b = (2, 4)
    black_stones.add(move1_b)
    analyze_and_print(1, move1_b, black_stones, white_stones)

    # 2. W:(1,5) - Forced response
    move2_w = (1, 5)
    white_stones.add(move2_w)
    analyze_and_print(2, move2_w, black_stones, white_stones)

    # 3. B:(1,6)
    move3_b = (1, 6)
    black_stones.add(move3_b)
    analyze_and_print(3, move3_b, black_stones, white_stones)

    # 4. W:(1,3) - Forced response
    move4_w = (1, 3)
    white_stones.add(move4_w)
    analyze_and_print(4, move4_w, black_stones, white_stones)

    # 5. B:(2,3) - The key move (Tesuji)
    move5_b = (2, 3)
    black_stones.add(move5_b)
    analyze_and_print(5, move5_b, black_stones, white_stones)
    print("This is a DOUBLE ATARI. White can only save one group.")

    # 6. W:(3,2) - White chooses to save the bottom group
    move6_w = (3, 2)
    white_stones.add(move6_w)
    analyze_and_print(6, move6_w, black_stones, white_stones)
    print("White saves the bottom group, but the top group will be captured.")
    
    # 7. B:(1,2) - Black captures the top group
    move7_b = (1, 2)
    captured_group = find_group(1, 3, white_stones)
    white_stones -= captured_group
    black_stones.add(move7_b)
    print(f"--- Turn 7: Black plays at {move7_b} ---")
    print(f"Black captures the White group at: {', '.join(map(str, sorted(list(captured_group))))}")
    
    # Final analysis of the remaining white stones
    analyze_and_print(7, "Post-Capture", black_stones, white_stones)
    print("The single remaining White group is surrounded with only 2 liberties and will be captured.")

if __name__ == '__main__':
    solve_go_problem()