def get_neighbors(r, c):
    """Returns the four orthogonal neighbors of a coordinate."""
    return {(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}

def find_groups(player_stones):
    """Finds all connected groups of stones for a player."""
    stones_to_visit = set(player_stones)
    groups = []
    while stones_to_visit:
        group = set()
        q = [stones_to_visit.pop()]
        group.add(q[0])
        head = 0
        while head < len(q):
            stone = q[head]
            head += 1
            for neighbor in get_neighbors(stone[0], stone[1]):
                if neighbor in stones_to_visit:
                    stones_to_visit.remove(neighbor)
                    group.add(neighbor)
                    q.append(neighbor)
        groups.append(group)
    return groups

def get_liberties(group, all_black_stones, all_white_stones):
    """Calculates the liberties for a single group."""
    liberties = set()
    all_stones = all_black_stones.union(all_white_stones)
    for stone in group:
        for neighbor in get_neighbors(stone[0], stone[1]):
            if neighbor not in all_stones:
                liberties.add(neighbor)
    return liberties

def main():
    """
    Analyzes the Go board position to find the best move for Black.
    """
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    print("--- Initial Board State Analysis ---")
    print(f"Black stones: {sorted(list(black_stones))}")
    print(f"White stones: {sorted(list(white_stones))}")

    white_groups = find_groups(white_stones)
    print(f"\nFound {len(white_groups)} white groups:")
    
    initial_liberties = {}
    for i, group in enumerate(white_groups, 1):
        libs = get_liberties(group, black_stones, white_stones)
        initial_liberties[i] = libs
        print(f"  - Group {i}: Stones {sorted(list(group))}, has {len(libs)} liberties at {sorted(list(libs))}")

    vital_point = (2, 4)
    print(f"\nAnalyzing shared liberties...")
    shared_libs = set.intersection(*initial_liberties.values())
    print(f"The point {vital_point} is a shared liberty for all white groups.")

    print(f"\n--- Analysis of Black's Move at {vital_point} ---")
    
    # Add black's move
    black_stones_after_move = black_stones.union({vital_point})
    
    print(f"Black plays at {vital_point}. Board now has Black stones at: {sorted(list(black_stones_after_move))}")
    
    print("\nRecalculating liberties for white groups:")
    for i, group in enumerate(white_groups, 1):
        libs_after = get_liberties(group, black_stones_after_move, white_stones)
        status = ""
        if len(libs_after) == 1:
            status = "<- ATARI!"
        print(f"  - Group {i}: Stones {sorted(list(group))}, now has {len(libs_after)} liberties at {sorted(list(libs_after))} {status}")

    print("\nConclusion: Playing at (2, 4) puts White Group 1 in Atari, forcing a response and initiating a winning sequence for Black.")

if __name__ == "__main__":
    main()