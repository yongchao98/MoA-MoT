def get_sharps(note_name):
    """
    Calculates the number of sharps in a major key signature based on the note name,
    following the problem's specific rules.
    """
    if not isinstance(note_name, str) or not note_name:
        raise ValueError("Note name must be a non-empty string")

    # Rule: F major (1 flat) is treated as having 11 sharps.
    if note_name == "F":
        return 11

    # Base values for natural notes (except F) according to the Circle of Fifths.
    base_sharps_map = {
        'C': 0, 'G': 1, 'D': 2, 'A': 3, 'E': 4, 'B': 5
    }

    letter = note_name[0]
    num_sharps = note_name.count('#')

    if letter not in base_sharps_map and letter != 'F':
        raise ValueError(f"Invalid note letter: {letter}")

    # Rule: Special handling for F#, F##, etc.
    if letter == 'F':
        # F# major has 6 sharps.
        # Each subsequent sharp adds 7.
        return 6 + (num_sharps - 1) * 7

    # Rule: For all other notes, adding a sharp to the tonic adds 7 sharps to the signature.
    base_val = base_sharps_map[letter]
    return base_val + num_sharps * 7

def main():
    """
    Derives the formula by calculating sums for n=0, 1, 2 and finding the pattern.
    """
    initial_notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]

    sums = []
    print("Step 1: Calculate the sum of sharps for n=0, 1, and 2.")
    # Calculate for n=0, n=1, and n=2 to establish a pattern
    for n in range(3):
        current_notes = [note + '#' * n for note in initial_notes]
        total_sharps = sum(get_sharps(note) for note in current_notes)
        sums.append(total_sharps)
        print(f"For n={n}, the total sum of sharps is: {total_sharps}")

    sum_n0, sum_n1, sum_n2 = sums

    print("\nStep 2: Analyze the pattern in the sums.")
    # The change from n=0 to n=1 is anomalous due to the F -> F# transition
    diff_0_1 = sum_n1 - sum_n0
    print(f"The difference between Sum(1) and Sum(0) is {sum_n1} - {sum_n0} = {diff_0_1}.")

    # The change for n>0 is constant because the note 'F' no longer appears in the tonic list.
    diff_1_2 = sum_n2 - sum_n1
    print(f"The difference between Sum(2) and Sum(1) is {sum_n2} - {sum_n1} = {diff_1_2}.")
    print(f"For n > 0, the sum increases by a constant {diff_1_2} for each step.")
    
    print("\nStep 3: Derive the formula for n > 0.")
    print("The formula is an arithmetic progression for n > 0 of the form: Sum(n) = Sum(1) + (n-1) * d")
    
    # The common difference 'd' for n>0 is diff_1_2
    d = diff_1_2
    # The starting term for the progression (n>0) is Sum(1)
    s1 = sum_n1
    
    # The formula is Sum(n) = s1 + (n-1)*d = s1 + d*n - d = (s1 - d) + d*n
    constant_term = s1 - d
    n_coefficient = d
    
    print(f"Sum(n) = {s1} + (n-1) * {d}")
    print(f"Sum(n) = {s1} + {d}*n - {d}")
    print(f"Simplified formula for n > 0: Sum(n) = {constant_term} + {n_coefficient}*n")
    
    print("\nFinal Answer:")
    print("The simplified formula for the sum of the number of sharps for n > 0 is:")
    print(f"{constant_term} + {n_coefficient}n")


if __name__ == "__main__":
    main()