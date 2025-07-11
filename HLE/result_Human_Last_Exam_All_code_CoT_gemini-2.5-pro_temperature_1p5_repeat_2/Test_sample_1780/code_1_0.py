import collections

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters L, N, S, W
    based on a specific connection rule for their tips.
    """
    # Step 1: Define letters and their tip positions (Left, Right)
    letter_tips = {
        'L': ('top', 'bottom'),
        'N': ('bottom', 'top'),
        'S': ('top', 'bottom'),
        'W': ('top', 'top')
    }
    all_letters = list(letter_tips.keys())

    # Step 2: Build a graph based on the connection rule.
    # An edge from u -> v exists if right_tip(u) is different from left_tip(v).
    graph = collections.defaultdict(list)
    for u in all_letters:
        for v in all_letters:
            if u == v:
                continue
            
            u_right_tip = letter_tips[u][1]
            v_left_tip = letter_tips[v][0]
            
            # The connection rule derived from the example (S->W works, W->S doesn't)
            if u_right_tip != v_left_tip:
                graph[u].append(v)

    # Step 3: Find all Hamiltonian paths (arrangements of all 4 letters) using DFS.
    found_paths = []

    def find_all_paths(current_path):
        """Recursively search for paths that visit every letter exactly once."""
        if len(current_path) == len(all_letters):
            found_paths.append(current_path)
            return

        last_letter_in_path = current_path[-1]
        for neighbor in graph[last_letter_in_path]:
            if neighbor not in current_path:
                find_all_paths(current_path + [neighbor])

    for start_letter in all_letters:
        find_all_paths([start_letter])

    # Step 4: Output the results and the final count.
    print("The valid arrangements found are:")
    for i, path in enumerate(found_paths):
        print(f"  {i+1}. {''.join(path)}")
        
    print("\nThe final calculation is based on the count of each valid path:")
    
    if not found_paths:
        total = 0
        print("0 = 0")
    else:
        # For the prompt "output each number in the final equation", we represent
        # the total as a sum of 1 for each path found.
        sum_components = ["1"] * len(found_paths)
        total = len(found_paths)
        print(f"{' + '.join(sum_components)} = {total}")

    print(f"\nThe total number of ways is {total}.")

solve_letter_arrangement()