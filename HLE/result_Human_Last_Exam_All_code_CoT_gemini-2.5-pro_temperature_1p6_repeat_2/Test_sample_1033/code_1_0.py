import string

def solve():
    """
    Calculates the next three terms in the sequence based on a derived formula.
    """
    # The last index in the given sequence is 83. We need to find terms for indices 84, 85, 86.
    start_index = 84
    num_terms_to_predict = 3

    # As per the plan, L1 for the next terms will be 'O'.
    l1 = 'O'
    v1 = ord(l1) - ord('A')

    # As per the plan, L2 will start from 'A' and increment.
    l2_start_char = 'A'
    
    print("The formula is: L3 = (16 * value(L2) + 14 * index + 21) % 26\n")

    results = []
    for i in range(num_terms_to_predict):
        current_index = start_index + i
        
        # Determine L2 and its value v2
        l2 = chr(ord(l2_start_char) + i)
        v2 = ord(l2) - ord('A')
        
        # Apply the formula to find v3
        v3 = (16 * v2 + 14 * current_index + 21) % 26
        l3 = chr(ord('A') + v3)
        
        # The complete triplet
        triplet = f"{l1}{l2}{l3}"
        results.append(triplet)
        
        # Print the equation for the current term
        print(f"For index {current_index} (L1='{l1}', L2='{l2}'):")
        print(f"(16 * {v2} + 14 * {current_index} + 21) % 26 = {v3}")
        print(f"Resulting triplet: {triplet}\n")

    final_answer = " ".join(results)
    print(f"The next three capital letters in the sequence are: {final_answer}")


solve()
<<<OAB OBF OCJ>>>