import random

def alice_strategy_A(sequence, box_k_index, modulus):
    """
    Simulates Alice's strategy for a sequence in case (A).

    Args:
        sequence (list): The sequence of numbers in the boxes (must have finite support).
        box_k_index (int): The index of the box Alice leaves closed (1-based).
        modulus (int): The modulus for the arithmetic, e.g., 10 for 1/10 prob.
    """
    print(f"--- Alice's Strategy Simulation ---")
    print(f"The true sequence is (partially): {sequence[:box_k_index+10]}...")
    
    # 1. Alice leaves box k closed and opens all others.
    closed_box_value = sequence[box_k_index - 1]
    print(f"Alice leaves box #{box_k_index} closed. The true value is {closed_box_value}.")

    # 2. She computes the sum of all boxes she opened.
    partial_sum = sum(sequence) - closed_box_value
    print(f"Alice opens all other boxes and calculates their sum: S' = {partial_sum}")

    # 3. Alice bets on the equivalence class of the total sum.
    # She picks a class j* at random from {0, 1, ..., modulus-1}.
    j_star = random.randint(0, modulus - 1)
    print(f"Alice randomly bets that the total sum modulo {modulus} will be: j* = {j_star}")

    # 4. Alice formulates her guess 'g' based on her bet.
    # The equation she solves is: (g + S') mod modulus = j*
    # g mod modulus = (j* - S') mod modulus
    # She makes the simplest guess, the smallest non-negative integer.
    guess = (j_star - partial_sum) % modulus
    
    print("\n--- The Final Equation ---")
    print(f"Alice's guess 'g' is based on the equation:")
    print(f"(g + {partial_sum}) mod {modulus} = {j_star}")
    print(f"Solving for g: g = ({j_star} - {partial_sum}) mod {modulus}")
    print(f"Alice's final guess for the number in box #{box_k_index} is: {guess}")

    # 5. Check if she succeeds.
    # Her success depends on her bet j* being correct (and s(k) being small).
    total_sum = sum(sequence)
    j_true = total_sum % modulus
    
    print("\n--- Result ---")
    print(f"The true total sum is {total_sum}, so the true class is j_true = {j_true}.")
    if guess == closed_box_value:
        print(f"Success! Alice's guess ({guess}) matches the true value ({closed_box_value}).")
    else:
        print(f"Failure. Alice's guess ({guess}) does not match the true value ({closed_box_value}).")
    
    print(f"Success was only guaranteed if her random bet j*={j_star} matched the true class j_true={j_true} (and other assumptions held).")
    print(f"The probability of this is 1/{modulus}.")


# Example execution for Case (A)
# Create a sequence with finite support (all zeros after some index).
s_A = [5, 8, 0, 2, 4, 1, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0] 
# To get P >= 9/10, we would need a modulus <= 10/9, which is impossible.
# We demonstrate with modulus 10, which gives P=1/10, showing a valid strategy exists.
alice_strategy_A(sequence=s_A, box_k_index=5, modulus=10)
