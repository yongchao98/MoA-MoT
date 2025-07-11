import math

def calculate_probability():
    """
    Calculates the probability that for each individual, there exists a type of item 
    for which they hold strictly more copies than any other individual.
    """

    # Step 1: Calculate S, the total number of ways to distribute the items.
    # S = 25! / (5! * 5! * 5! * 5! * 5!)
    f5 = math.factorial(5)
    S = math.factorial(25) // (f5**5)

    # Step 2: Calculate F, the number of favorable distributions.
    # F = 5! * F_1, where F_1 is the number of ways for a fixed person-type dominance.
    # We calculate F_1 by casework on the number of "P4" columns, i.e.,
    # columns of the form (5, 0, 0, 0, 0).

    # Case 0: 0 P4-columns. This implies the count matrix has diagonals of 2.
    # There are 6 such matrices (related to 3-regular digraphs on 5 vertices).
    # For each, W(C) = (5! / (2!*1!*1!*1!*0!))^5 = 60^5
    f1_c0 = 6 * (f5 // 2)**5

    # Case 1: 1 P4-column. 5 choices for the column.
    # The subproblem is a 4x4 matrix with row/col sums of 5.
    # The unique solution for the submatrix has diagonals of 2.
    # W(C) = (5!/5!) * (5!/2!)^4 = 1 * 60^4
    f1_c1 = 5 * (f5 // 2)**4

    # Case 2: 2 P4-columns. C(5,2)=10 choices.
    # Subproblem is 3x3 with row/col sums of 5. Unique solution has diagonals of 3.
    # W(C) = (5!/5!)^2 * (5!/3!)^3 = 1^2 * 20^3
    f1_c2 = 10 * (f5 // 6)**3

    # Case 3: 3 P4-columns. C(5,3)=10 choices.
    # Subproblem is 2x2 with row/col sums of 5. Two solutions without P4-like columns.
    # Sol 1: Diagonals are 4,1. W(C) = (5!/5!)^3 * (5!/4!)^2 = 5^2 = 25
    # Sol 2: Diagonals are 3,2. W(C) = (5!/5!)^3 * (5!/(3!2!))^2 = 10^2 = 100
    n_c3_sol1 = (f5 // math.factorial(4))**2
    n_c3_sol2 = (f5 // (math.factorial(3) * math.factorial(2)))**2
    f1_c3 = 10 * (n_c3_sol1 + n_c3_sol2)

    # Case 4: 4 P4-columns. C(5,4)=5 choices.
    # Subproblem is 1x1 with sum 5. Unique solution [5]. W(C) = 1.
    f1_c4 = 5 * 1

    # Case 5: 5 P4-columns. C(5,5)=1 choice.
    # The matrix is 5*Identity. W(C) = 1.
    f1_c5 = 1 * 1

    # Summing all cases for F_1
    F_1 = f1_c0 + f1_c1 + f1_c2 + f1_c3 + f1_c4 + f1_c5
    
    # Total favorable ways F
    F = f5 * F_1
    
    # Step 3: Calculate the probability P = F / S
    P = F / S

    print(f"Total number of distributions (S):")
    print(S)
    print("\nNumber of favorable distributions (F):")
    print(F)
    print("\nProbability P = F / S:")
    print(f"{F} / {S} = {P}")
    
    return P

# Execute the calculation and store the final probability.
final_probability = calculate_probability()

# The final answer in the specified format.
# The calculation shows P is approximately 0.00091065.
# Let's print the result in the required format.
print(f"\n<<< {final_probability} >>>")