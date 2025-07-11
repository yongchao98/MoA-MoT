import math

def calculate_and_print_solution():
    """
    Calculates the imaginary part of the sum of the integrals and prints the steps.
    """
    # Step 1: Define the values of the residues
    res_1_5 = 3 * math.sqrt(math.pi) / 4
    res_neg_1 = -2 / 5
    res_neg_2 = 2 / 7
    res_neg_3 = -1 / 9

    # Step 2: Sum the winding numbers for each pole
    w_total_1_5 = 1 + 0
    w_total_neg_1 = 1 + 1
    w_total_neg_2 = -1 + 1
    w_total_neg_3 = 1 + 0
    
    # Step 3: Calculate the final sum
    final_sum_in_brackets = (w_total_1_5 * res_1_5 + 
                             w_total_neg_1 * res_neg_1 +
                             w_total_neg_2 * res_neg_2 + 
                             w_total_neg_3 * res_neg_3)
                             
    final_result = 2 * math.pi * final_sum_in_brackets

    # Step 4: Print the full explanation and calculation
    print("The imaginary part of the sum of the integrals is S.")
    print("S = 2 * pi * Sum[ (W(C1, zk) + W(C2, zk)) * Res(f, zk) ]")
    print("\nPlugging in the winding numbers and residues:")
    print("S = 2 * pi * [ (1+0) * (3*sqrt(pi)/4) + (1+1) * (-2/5) + (-1+1) * (2/7) + (1+0) * (-1/9) ]")
    print(f"S = 2 * pi * [ {w_total_1_5} * ({res_1_5:.6f}) + {w_total_neg_1} * ({res_neg_1:.6f}) + {w_total_neg_2} * ({res_neg_2:.6f}) + {w_total_neg_3} * ({res_neg_3:.6f}) ]")
    print("The contribution from the pole at z = -2 cancels out.")
    print("S = 2 * pi * [ 3*sqrt(pi)/4 - 4/5 - 1/9 ]")
    print(f"\nNumerically:")
    print(f"S = 2 * {math.pi:.6f} * [ {res_1_5:.6f} - {4/5:.6f} - {abs(res_neg_3):.6f} ]")
    print(f"S = 2 * {math.pi:.6f} * [ {final_sum_in_brackets:.6f} ]")
    print(f"S = {final_result:.6f}")

calculate_and_print_solution()