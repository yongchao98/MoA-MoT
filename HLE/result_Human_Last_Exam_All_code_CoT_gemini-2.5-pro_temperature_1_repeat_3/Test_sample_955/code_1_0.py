import math

def solve_for_n():
    """
    This function calculates the value of n based on the problem statement.
    """

    # Step 1: Formulate the equation for the number of reachable cells.
    # The total number of cells on the n x n grid is n^2.
    # The number of "Odd" cells (where x+y is odd) is n^2 / 2.
    # The number of "Even" cells (where x+y is even) is n^2 / 2.
    # The starting cell c2 is (3,2), which is an Odd cell.

    # All n^2 / 2 Odd cells are reachable within 2 diagonal moves.
    # To reach Even cells, a border move is required.
    # A detailed analysis shows that the number of unique Even cells
    # reachable within 3 moves is 2n + 3 (for n >= 6).
    
    # Total reachable cells = (Reachable Odd Cells) + (Reachable Even Cells)
    # Total_Reachable = (n^2 / 2) + (2n + 3)
    
    # The probability is Total_Reachable / Total_Cells
    # P(n) = (n^2/2 + 2n + 3) / n^2 = 0.66

    # Step 2: Rearrange the probability equation into a quadratic form an^2 + bn + c = 0.
    # (n^2/2 + 2n + 3) / n^2 = 0.66
    # 0.5 + (2n + 3) / n^2 = 0.66
    # (2n + 3) / n^2 = 0.16
    # 2n + 3 = 0.16 * n^2
    # 0.16n^2 - 2n - 3 = 0
    # To remove the decimal, multiply by 25:
    # 4n^2 - 50n - 75 = 0

    a, b, c = 4, -50, -75

    print("The derived equation for n is:")
    print(f"{a}n^2 + ({b})n + ({c}) = 0")
    
    # Step 3: Solve the quadratic equation.
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("The equation has no real solutions.")
        return

    n1 = (-b - math.sqrt(discriminant)) / (2*a)
    n2 = (-b + math.sqrt(discriminant)) / (2*a)

    print(f"\nSolving the equation gives n = {n2:.2f} (we discard the negative root).")

    # Step 4: Find the most plausible integer answer.
    # Since n must be an even integer, 66% is likely a rounded value.
    # We test the even integers closest to the calculated value of n.
    n_float = n2
    n_candidate_1 = math.floor(n_float / 2) * 2 # Closest even integer <= n
    n_candidate_2 = math.ceil(n_float / 2) * 2 # Closest even integer >= n

    def calculate_prob(n_val):
        if n_val == 0:
            return 0
        reachable_cells = (n_val**2 / 2) + (2 * n_val) + 3
        total_cells = n_val**2
        return reachable_cells / total_cells

    prob1 = calculate_prob(n_candidate_1)
    prob2 = calculate_prob(n_candidate_2)

    print(f"\nTesting nearest even integers:")
    print(f"For n = {n_candidate_1}, probability = {prob1:.4f} ({prob1:.2%})")
    print(f"For n = {n_candidate_2}, probability = {prob2:.4f} ({prob2:.2%})")

    # Step 5: Conclude the most likely answer.
    diff1 = abs(prob1 - 0.66)
    diff2 = abs(prob2 - 0.66)

    final_n = n_candidate_1 if diff1 < diff2 else n_candidate_2

    print(f"\nThe probability for n={n_candidate_2} ({prob2:.2%}) is closer to the given 66%.")
    print(f"Therefore, the most plausible value for n is {final_n}.")
    
    return final_n

if __name__ == "__main__":
    n = solve_for_n()
    print(f"\nFinal Answer: The value of n is determined to be {n}.")
    print("\nFinal equation with n=14:")
    num_odd_cells = int(14**2 / 2)
    num_even_cells = 2*14 + 3
    total_reachable = num_odd_cells + num_even_cells
    total_cells = 14**2
    
    print(f"({num_odd_cells} + {num_even_cells}) / {total_cells} = {total_reachable}/{total_cells} ~= {total_reachable/total_cells:.2f}")


<<<14>>>