import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an isosceles right triangle.
    """
    # Explanation of the principle
    print("This solution uses a principle from geometric probability:")
    print("P(escape through a side) = (Length of that side) / (Total Perimeter)")
    print("-" * 50)

    # Step 1: Define the geometry. The choice of L is arbitrary as it will cancel out. We can use L=1.
    L = 1.0
    leg_length = L
    hypotenuse_length = L * math.sqrt(2)

    print("Step 1: Define side lengths for an isosceles right triangle (let leg length be L):")
    print(f"Length of each leg = L")
    print(f"Length of hypotenuse = L * sqrt(2) ≈ {hypotenuse_length:.4f} * L\n")

    # Step 2: Calculate the perimeter
    perimeter = leg_length + leg_length + hypotenuse_length
    print("Step 2: Calculate the total perimeter:")
    print(f"Perimeter = L + L + L * sqrt(2) = L * (2 + sqrt(2)) ≈ {perimeter:.4f} * L\n")

    # Step 3: Calculate the individual probabilities
    # Probability of escaping through the hypotenuse, P(H)
    # The formula is (L * sqrt(2)) / (L * (2 + sqrt(2))) which simplifies to sqrt(2) - 1
    p_h = math.sqrt(2) - 1

    # Probability of escaping through the two legs combined, P(L)
    # The formula is (2 * L) / (L * (2 + sqrt(2))) which simplifies to 2 - sqrt(2)
    p_legs = 2 - math.sqrt(2)

    print("Step 3: Calculate the probabilities:")
    print("P(H) = (Hypotenuse Length) / Perimeter = (L * sqrt(2)) / (L * (2 + sqrt(2)))")
    print(f"     = sqrt(2) - 1 ≈ {p_h:.4f}")
    print("\nP(L) = (Combined Leg Length) / Perimeter = (2 * L) / (L * (2 + sqrt(2)))")
    print(f"     = 2 - sqrt(2) ≈ {p_legs:.4f}\n")
    
    # Step 4: Calculate the final difference
    difference = p_h - p_legs
    
    print("Step 4: Calculate the difference P(H) - P(L):")
    print("\nThe final equation is:")
    
    # Printing each number in the final equation as requested
    # The equation is (sqrt(2) - 1) - (2 - sqrt(2)) = 2*sqrt(2) - 3
    sqrt2_val = math.sqrt(2)
    print(f"({sqrt2_val:.4f} - 1.0000) - (2.0000 - {sqrt2_val:.4f}) = {difference:.4f}")
    
    print("\nSymbolically, the equation simplifies to:")
    print(f"(sqrt(2) - 1) - (2 - sqrt(2)) = 2*sqrt(2) - 3")
    
    print("\nFinal Answer:")
    print(f"{difference}")

solve_escape_probability()