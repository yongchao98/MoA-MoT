def explain_and_demonstrate_mf_dimensions():
    """
    This function explains and demonstrates the dimensional structure of
    Type-1, Type-2, and Type-3 fuzzy membership functions (MFs).
    """

    print("Analyzing the dimensional structure of fuzzy membership functions:")

    # --- Type-1 MF ---
    # Domain has 1 variable: x
    # We can represent it as: μ(x)
    def type_1_mf(x):
        # A simple example calculation
        return max(0, 1 - abs(x - 5) / 5)

    x_val = 7
    print(f"\n1. A Type-1 MF is defined on a 1-variable domain.")
    print(f"   - Equation form: μ(x)")
    print(f"   - Example call: μ(x={x_val})")
    # Show the "equation" with the number plugged in
    print(f"   - The function μ({x_val}) maps a single input to a single output.")


    # --- Type-2 MF ---
    # Domain has 2 variables: x, u
    # We can represent it as: μ(x, u)
    def type_2_mf(x, u):
        # Conceptual: returns a secondary grade. Let's imagine a simple form.
        return max(0, 1 - abs(x-5)/5 - abs(u-0.8)/0.2)

    u_val = 0.9
    print(f"\n2. A Type-2 MF is defined on a 2-variable domain.")
    print(f"   - Equation form: μ(x, u)")
    print(f"   - 'u' is the primary membership variable.")
    print(f"   - Example call: μ(x={x_val}, u={u_val})")
    # Show the "equation" with the numbers plugged in
    print(f"   - The function μ({x_val}, {u_val}) maps two inputs to an output (a secondary grade).")

    # --- Type-3 MF ---
    # Domain has 3 variables: x, u, v
    # We can represent it as: μ(x, u, v)
    def type_3_mf(x, u, v):
        # Conceptual: returns a tertiary grade.
        return max(0, 1 - abs(x-5)/5 - abs(u-0.8)/0.2 - abs(v-0.7)/0.3)

    v_val = 0.6
    print(f"\n3. A Type-3 MF is defined on a 3-variable domain.")
    print(f"   - Equation form: μ(x, u, v)")
    print(f"   - 'v' is the secondary membership variable.")
    print(f"   - Example call: μ(x={x_val}, u={u_val}, v={v_val})")
    # Show the "equation" with the numbers plugged in
    print(f"   - The function μ({x_val}, {u_val}, {v_val}) maps three inputs to an output (a tertiary grade).")


    print("\n" + "="*50)
    print("CONCLUSION:")
    print("The fundamental difference in dimensional structure from Type-2 to Type-3")
    print("is that the function's domain is expanded from two variables to three variables.")
    print("="*50)

if __name__ == "__main__":
    explain_and_demonstrate_mf_dimensions()