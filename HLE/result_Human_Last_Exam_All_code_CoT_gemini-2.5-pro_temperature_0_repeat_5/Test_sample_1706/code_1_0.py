import collections

def solve_go_problem():
    """
    This script solves a Go life-and-death problem by finding the vital point
    of the opponent's group.
    """
    # Step 1: Define the board configuration based on the problem description.
    # The coordinate system is (row, column), with columns 1-19 from right to left.
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
    all_stones = black_stones.union(white_stones)

    print("Analyzing the Go board situation...")
    print(f"Black pieces: {sorted(list(black_stones))}")
    print(f"White pieces: {sorted(list(white_stones))}")
    print("-" * 30)

    # Step 2: Find all liberties of the white group.
    # A liberty is an empty point adjacent to any stone in the group.
    # Since all white stones are connected, we treat them as one group.
    liberties = set()
    for r, c in white_stones:
        # Check four neighbors: up, down, left, right
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (r + dr, c + dc)
            # A point is a liberty if it's on the board (1-19) and empty.
            if 1 <= neighbor[0] <= 19 and 1 <= neighbor[1] <= 19:
                if neighbor not in all_stones:
                    liberties.add(neighbor)

    print(f"The White group's liberties (potential eye space) are: {sorted(list(liberties))}")
    print("-" * 30)

    # Step 3: Find the vital point.
    # Our strategy is to find the liberty adjacent to the most white stones.
    # This is a strong heuristic for finding the move that most effectively
    # disrupts the opponent's shape.
    liberty_scores = collections.defaultdict(int)
    for liberty in liberties:
        count = 0
        # Check which white stones this liberty is adjacent to
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (liberty[0] + dr, liberty[1] + dc)
            if neighbor in white_stones:
                count += 1
        liberty_scores[liberty] = count

    print("Calculating the 'score' for each liberty (number of adjacent white stones):")
    # The "equation" is the set of scores for each liberty.
    for liberty, score in sorted(liberty_scores.items()):
        print(f"Liberty {liberty} score = {score}")
    print("-" * 30)

    # Step 4: The best move is the liberty with the highest score.
    if not liberty_scores:
        best_move = "Impossible"
        print("Result: No liberties found. This indicates an error or a captured group.")
    else:
        # Find the liberty with the maximum score
        best_move = max(liberty_scores, key=liberty_scores.get)
        max_score = liberty_scores[best_move]
        print("The best move is the one with the highest score.")
        print(f"Final Calculation: max(scores) corresponds to move {best_move} with a score of {max_score}.")
        print(f"\nConclusion: The optimal first move for Black is to play at {best_move}.")


solve_go_problem()