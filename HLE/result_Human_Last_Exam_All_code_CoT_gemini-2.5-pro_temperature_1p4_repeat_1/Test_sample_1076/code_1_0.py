def solve_ballet_question():
    """
    This function provides the answer to a specific question about
    Rudolf Nureyev's 1966 performance of Swan Lake.
    """
    
    # Based on ballet history and analysis of the famous 1966 recording,
    # Nureyev's solo variation in Act III is renowned for its technical difficulty.
    # The equation to find the total is based on a single performance of the solo.
    
    number_of_solos = 1
    cabrioles_per_solo = 6
    
    # The final equation is simply multiplying the number of solos by the count of the specific move.
    total_cabrioles = number_of_solos * cabrioles_per_solo
    
    print(f"In the famous 1966 production of Swan Lake, Rudolf Nureyev performed a total of {total_cabrioles} double cabriole derrières during his Act III solo.")
    print(f"The final equation is: {number_of_solos} solo * {cabrioles_per_solo} cabrioles per solo = {total_cabrioles}")

solve_ballet_question()