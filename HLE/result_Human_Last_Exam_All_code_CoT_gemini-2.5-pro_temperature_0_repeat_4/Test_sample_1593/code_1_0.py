import math

def solve_complexity():
    """
    This function analyzes the query complexity for two regimes and prints the result.
    """
    
    print("### Query Complexity Analysis ###\n")
    
    # --- Step 1 & 2: General Complexity Analysis ---
    print("Step 1: Two main algorithms are considered:")
    print("  a) Comparison Sort (e.g., MergeSort): Uses C-queries on full strings.")
    print("     - Complexity: Theta(N * log N)")
    print("  b) Radix Sort: Uses H-queries on chunks and C-queries to sort unique chunks.")
    print("     - Optimized Complexity: Theta(N * L / log N)")
    print("\nThe overall complexity is the minimum of the two: min(Theta(N*log N), Theta(N*L/log N)).\n")

    # --- Step 3: Regime 1 Analysis ---
    print("--- Regime 1: N = 2^sqrt(L) ---")
    print("From the relation, we derive L in terms of N: L = (log N)^2.")
    print("Substituting L into the radix sort complexity:")
    print("  Q_radix = Theta(N * (log N)^2 / log N) = Theta(N * log N).")
    print("Both algorithms yield Theta(N * log N) complexity.")
    
    # --- Step 4: Regime 2 Analysis ---
    print("\n--- Regime 2: N = 2^((log L)^2) ---")
    print("From the relation, we derive L in terms of N: L = 2^sqrt(log N).")
    print("Substituting L into the radix sort complexity:")
    print("  Q_radix = Theta(N * 2^sqrt(log N) / log N).")
    print("Comparing Q_merge = Theta(N*log N) with Q_radix:")
    print("  We compare (log N) with (2^sqrt(log N) / log N).")
    print("  Asymptotically, 2^sqrt(log N) grows faster than (log N)^2.")
    print("  Therefore, Q_radix is larger, and the minimum complexity is Theta(N * log N) from comparison sort.")

    # --- Step 5: Conversion to (a,b,c) notation ---
    print("\n--- Conversion to (a,b,c) Notation ---")
    print("For both regimes, the complexity is Theta(N * log N).")
    print("We need to find a, b, c such that: sqrt(N^a * (log N)^b * (log log N)^c) = N * log N")
    print("Squaring both sides: N^a * (log N)^b * (log log N)^c = N^2 * (log N)^2")
    
    a = 2
    b = 2
    c = 0
    
    print(f"By comparing exponents, we get:")
    print(f"  a = {a}")
    print(f"  b = {b}")
    print(f"  c = {c}")
    
    regime1_complexity = (a, b, c)
    regime2_complexity = (a, b, c)
    
    final_answer_string = f"{regime1_complexity},{regime2_complexity}".replace(" ", "")
    
    print("\n### Final Answer ###")
    print(f"The final answer in the format (a,b,c),(d,e,f) is: {final_answer_string}")

solve_complexity()