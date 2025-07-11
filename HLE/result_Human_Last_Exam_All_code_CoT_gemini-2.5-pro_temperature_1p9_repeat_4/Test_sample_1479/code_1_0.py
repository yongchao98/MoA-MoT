def explain_fuzzy_mf_dimensions():
    """
    Illustrates the dimensional structure difference between Type-1, Type-2,
    and Type-3 fuzzy membership functions (MFs).
    """
    print("--- Analyzing Dimensional Structure of Fuzzy Membership Functions ---")

    # --- 1. Type-1 Membership Function ---
    # A Type-1 MF is 2-dimensional. It maps one input variable to one output (membership value).
    # It's a function mu = f(x).
    print("\n[1] Type-1 Membership Function (2-Dimensional):")
    print("    - Maps an input 'x' to a membership value 'mu'.")
    print("    - Structure: mu = f(x)")
    # Example: For a temperature 'x' = 25, the membership in 'Warm' might be mu = 0.8.
    # The dimensions are (x, mu).
    x1 = 25
    mu1 = 0.8
    print(f"    - Example: f({x1}) = {mu1}")

    # --- 2. Type-2 Membership Function ---
    # A Type-2 MF is 3-dimensional. It models uncertainty about the membership value itself.
    # It maps an input 'x' AND a potential primary membership 'u' to a secondary membership value.
    # It's a function mu_secondary = f(x, u).
    print("\n[2] Type-2 Membership Function (3-Dimensional):")
    print("    - Maps an input 'x' and a primary membership 'u' to a secondary membership.")
    print("    - Structure: mu_secondary = f(x, u)")
    # Example: For temperature 'x' = 25, the belief that the primary membership 'u' = 0.8 is correct might be mu_secondary = 1.0.
    # The dimensions are (x, u, mu_secondary). This adds a 3rd dimension over Type-1.
    x2 = 25
    u2 = 0.8 # Primary membership value being evaluated
    mu_secondary = 1.0
    print(f"    - Example: f({x2}, {u2}) = {mu_secondary}")

    # --- 3. Type-3 Membership Function (The Core Difference) ---
    # A Type-3 MF is 4-dimensional. It introduces a TERTIARY membership function.
    # It models uncertainty about the secondary membership from the Type-2 set.
    # It maps 'x', primary membership 'u', AND secondary membership 'v' to a tertiary value.
    # It's a function mu_tertiary = f(x, u, v).
    print("\n[3] Type-3 Membership Function (4-Dimensional):")
    print("    - This is the key structural shift from Type-2.")
    print("    - It INTRODUCES A TERTIARY MEMBERSHIP FUNCTION to handle more complex uncertainty.")
    print("    - Structure: mu_tertiary = f(x, u, v)")
    # Example: For temp 'x'=25 and primary membership 'u'=0.8, the belief that the secondary membership 'v'=1.0 is correct might be mu_tertiary = 0.9.
    # The dimensions are (x, u, v, mu_tertiary).
    x3 = 25
    u3 = 0.8 # Primary membership
    v3 = 1.0 # Secondary membership being evaluated
    mu_tertiary = 0.9
    print(f"    - Example: f({x3}, {u3}, {v3}) = {mu_tertiary}")

    # --- Conclusion and Final Equation ---
    print("\n" + "="*60)
    print("Conclusion:")
    print("The fundamental structural difference from Type-2 to Type-3 is")
    print("the addition of another variable to the function, which requires the")
    print("introduction of a 'Tertiary Membership Function'. This adds a 4th dimension.")
    
    # Constraint: Print each number in a final equation.
    # The equation shows dimensionality increasing from Type-2 (3D) to Type-3 (4D).
    print("\nDimensionality Equation (Type-2 to Type-3):")
    type_2_dims = 3
    added_dims = 1
    type_3_dims = type_2_dims + added_dims
    
    print(f"{type_3_dims} = {type_2_dims} + {added_dims}")
    
    print("\nPrinting each number/symbol in the final equation:")
    print(type_3_dims)
    print("=")
    print(type_2_dims)
    print("+")
    print(added_dims)
    print("="*60)

# Run the explanation
explain_fuzzy_mf_dimensions()