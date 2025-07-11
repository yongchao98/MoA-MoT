import math

def solve_pm():
    """
    Calculates and prints the formula for P_m.
    This function doesn't take an argument as it prints the derived symbolic formula.
    """
    print("To find the probability P_m, we compute the number of favorable pairs (i, j) and divide by the total number of pairs.")
    print("\n1. Total Number of Pairs:")
    print("The total number of ways to choose two distinct indices from {1, 2, ..., 4m+2} is:")
    print("N_total = C(4m+2, 2) = (4m+2)(4m+1) / 2 = (2m+1)(4m+1)")

    print("\n2. Favorable Pairs (N_good):")
    print("The number of favorable pairs is determined by analyzing partitioning possibilities based on the parity of m.")
    print("For odd m, N_good = (m+1)(m+2) / 2")
    print("For even m, N_good = [(m+1)(m+2) / 2] + 1 = (m^2 + 3m + 4) / 2")
    print("\nA unified formula for N_good for any positive integer m can be expressed using (-1)^m:")
    print("N_good = (m^2 + 3m + 3 - (-1)^m) / 2")
    
    print("\n3. The Probability P_m:")
    print("P_m is the ratio N_good / N_total.")
    print("\nP_m = [ (m^2 + 3m + 3 - (-1)^m) / 2 ] / [ (2m+1)(4m+1) ]")
    print("\nThe final formula for P_m is:")
    print("      m^2 + 3m + 3 - (-1)^m")
    print("P_m = -------------------------")
    print("        2(2m+1)(4m+1)")

solve_pm()