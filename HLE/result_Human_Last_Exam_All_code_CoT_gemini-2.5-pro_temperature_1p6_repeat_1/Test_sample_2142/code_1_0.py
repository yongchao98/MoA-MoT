import requests
import sys

def solve_riddle():
    """
    Solves the Genesis P-Orridge riddle by finding words that match the given patterns and constraints.
    """
    # Step 1: Fetch a word list
    try:
        # A common word list available online
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        words = response.text.splitlines()
        print("Successfully downloaded word list.", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the word list. {e}", file=sys.stderr)
        print("Please check your internet connection or the URL.", file=sys.stderr)
        return

    # Step 2: Filter words for Y and Z patterns
    y_candidates = []
    z_candidates = []

    print("Filtering 6-letter and 9-letter words...", file=sys.stderr)
    for word in words:
        # Pattern for Y: 6 letters, word[1] == word[4]
        if len(word) == 6:
            # X1, X2, X3, X4, X5 must be unique letters for Y to be valid
            if word[1] == word[4] and len(set(word[0:4] + word[5:6])) == 5:
                y_candidates.append(word.upper())
        
        # Pattern for Z: 9 letters, word[1] == word[8]
        elif len(word) == 9:
            if word[1] == word[8]:
                z_candidates.append(word.upper())

    print(f"Found {len(y_candidates)} candidates for Y and {len(z_candidates)} candidates for Z.", file=sys.stderr)
    print("Searching for a matching pair...", file=sys.stderr)

    # Step 3: Check all pairs for matching constraints
    for y_word in y_candidates:
        for z_word in z_candidates:
            # Check cross-word constraints
            # Y[0](X1) == Z[2](X1)
            # Y[2](X3) == Z[0](X3)
            # Y[1](X2) == Z[7](X2)
            if y_word[0] == z_word[2] and y_word[2] == z_word[0] and y_word[1] == z_word[7]:
                
                # Assign all characters to X variables to check for 10 unique letters
                try:
                    x1 = y_word[0]
                    x2 = y_word[1]
                    x3 = y_word[2]
                    x4 = y_word[3]
                    x5 = y_word[5]
                    
                    x6 = z_word[1]
                    x7 = z_word[3]
                    x8 = z_word[4]
                    x9 = z_word[5]
                    x10 = z_word[6]

                    all_x = {x1, x2, x3, x4, x5, x6, x7, x8, x9, x10}

                    # Check for the 10 unique letter constraint
                    if len(all_x) == 10:
                        # A solution is found
                        final_y = y_word
                        final_z = z_word

                        # Print the mapping as required
                        print("\n--- Solution Found ---", file=sys.stderr)
                        print(f"X1 = {x1}, X2 = {x2}, X3 = {x3}, X4 = {x4}, X5 = {x5}", file=sys.stderr)
                        print(f"X6 = {x6}, X7 = {x7}, X8 = {x8}, X9 = {x9}, X10 = {x10}", file=sys.stderr)
                        
                        # Print the final Y and Z words
                        print(f"\nY = {final_y}")
                        print(f"Z = {final_z}")
                        
                        # Print the words separated by a space as the final answer
                        print(f"\nFinal combined answer:")
                        print(f"{final_y} {final_z}")
                        return

    print("No solution found in the provided dictionary.", file=sys.stderr)


if __name__ == "__main__":
    solve_riddle()