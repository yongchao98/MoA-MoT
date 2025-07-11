def display_infinite_product_solution():
    """
    This function derives and displays the closed-form expression for the infinite product
    Π_{n=3 to ∞} (1 - z^3/n^3).
    It prints each component and number of the final formula as requested.
    """
    print("The final expression for the infinite product is a fraction.")
    print("----------------------------------------------------------\n")

    # Numerator
    numerator = 1
    print(f"The numerator of the fraction is:")
    print(numerator)
    print("\n")

    # Denominator
    print("The denominator is a product of three main parts:")
    print("Part 1: The term corresponding to n=1, which is (1 - z^3 / d1)")
    d1 = 1**3
    print(f"The number d1 is 1^3, which is:")
    print(d1)
    print("\n")
    
    print("Part 2: The term corresponding to n=2, which is (1 - z^3 / d2)")
    d2 = 2**3
    print(f"The number d2 is 2^3, which is:")
    print(d2)
    print("\n")

    print("Part 3: The infinite product from n=1, expressed using the Gamma function Γ:")
    print("Γ(A - z) * Γ(A - z*ω) * Γ(A - z*ω^2)")
    A = 1
    print("The number A in the Gamma function arguments is:")
    print(A)
    print("\n")

    print("...where ω is the complex cube root of unity, ω = exp(B*π*i / C)")
    B = 2
    C = 3
    print("The number B in the definition of ω is:")
    print(B)
    print("The number C in the definition of ω is:")
    print(C)
    print("\n")

    # Final Assembled Equation
    print("----------------------------------------------------------")
    print("Assembling all parts, the final equation is:")
    print(f"\n      {numerator}")
    print("---------------------------------------------------------------------------------")
    print(f"(1 - z^3/{d1}) * (1 - z^3/{d2}) * Γ({A} - z) * Γ({A} - z*ω) * Γ({A} - z*ω^2)\n")
    print("Or, with simplified constants:")
    print(f"\n      {numerator}")
    print("-----------------------------------------------------------------")
    print(f"(1 - z^3) * (1 - z^3/8) * Γ(1 - z) * Γ(1 - z*ω) * Γ(1 - z*ω^2)")


display_infinite_product_solution()