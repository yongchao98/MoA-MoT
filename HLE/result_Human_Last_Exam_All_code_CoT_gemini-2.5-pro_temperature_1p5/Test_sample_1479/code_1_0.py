def explain_fuzzy_set_dimensions():
    """
    Explains the dimensional structure difference between Type-2 and Type-3 fuzzy sets.
    """
    print("Step 1: Understanding the baseline - The Type-1 Fuzzy Set")
    print("A Type-1 Fuzzy Membership Function (MF) is a 2D structure.")
    print(" - Variable 1: An element 'x' from the universe of discourse (the x-axis).")
    print(" - Variable 2: Its membership grade 'μ(x)' between 0 and 1 (the y-axis).")
    print("   Equation form: y = μ(x)")
    print("-" * 60)

    print("Step 2: Understanding the Type-2 Fuzzy Set")
    print("A Type-2 MF models uncertainty about the membership grade.")
    print("It expands the structure to three dimensions.")
    print(" - Variable 1: The element 'x' from the universe of discourse.")
    print(" - Variable 2: A primary membership variable 'u', representing possible membership grades for x.")
    print(" - Variable 3: A secondary membership grade 'μ(x, u)', which is the membership of 'u'.")
    print("This introduces a *Secondary Membership Function* to describe the uncertainty.")
    print("   The structure involves (x, u, μ(x, u)) making it 3D.")
    print("-" * 60)

    print("Step 3: Understanding the Type-3 Fuzzy Set")
    print("A Type-3 MF models uncertainty about the Type-2 secondary membership.")
    print("It expands the structure to four dimensions.")
    print(" - Variable 1: The element 'x' from the universe of discourse.")
    print(" - Variable 2: A primary membership variable 'u'.")
    print(" - Variable 3: A secondary membership variable 'v'.")
    print(" - Variable 4: A tertiary membership grade 'μ(x, u, v)', which is the membership of 'v'.")
    print("The fundamental change is the addition of a new layer of modeling, which requires a new function.")
    print("-" * 60)

    print("Conclusion: The Fundamental Difference")
    print("The transition from a Type-2 structure (3D) to a Type-3 structure (4D) is defined by the need for a new function to describe the newest layer of uncertainty.")
    print("This new function is called the 'Tertiary Membership Function'.")
    print("\nTherefore, the most accurate answer describing this fundamental structural change is:")
    print("F. Tertiary membership functions introduced")


if __name__ == '__main__':
    explain_fuzzy_set_dimensions()