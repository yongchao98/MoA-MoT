import sys

def explain_fgh_tripled_fixed_point_conditions():
    """
    This script outlines the conditions for the existence and uniqueness of an
    FGH-tripled fixed point based on the Banach Fixed-Point Theorem.
    """
    
    # Use Unicode for clarity if the terminal supports it
    use_unicode = sys.stdout.encoding.lower().startswith('utf')
    
    # Symbolic constants
    a1, a2, a3 = ("α₁", "α₂", "α₃") if use_unicode else ("a1", "a2", "a3")
    b1, b2, b3 = ("β₁", "β₂", "β₃") if use_unicode else ("b1", "b2", "b3")
    g1, g2, g3 = ("γ₁", "γ₂", "γ₃") if use_unicode else ("g1", "g2", "g3")
    
    print("--- Conditions for an FGH-Tripled Fixed Point ---\n")
    
    # Step 1: Define the mathematical space and functions
    print("Step 1: The Setting")
    print("1. Let (X, d_X), (Y, d_Y), and (Z, d_Z) be complete metric spaces.")
    print("   - A metric space is a set with a notion of distance (a metric 'd') between its elements.")
    print("   - 'Complete' means any Cauchy sequence of elements converges to a limit within the set.\n")
    print("2. Let F, G, and H be mappings defined as:")
    print("   - F: X * Y * Z -> X")
    print("   - G: Y * X * Y -> Y")
    print("   - H: Z * Y * X -> Z\n")
    
    # Step 2: Define what an FGH-tripled fixed point is
    print("Step 2: The Fixed-Point Equations")
    print("An FGH-tripled fixed point is a triplet of points (x, y, z) that simultaneously satisfies:")
    print("   F(x, y, z) = x")
    print("   G(y, x, y) = y")
    print("   H(z, y, x) = z\n")
    
    # Step 3: State the contraction conditions
    print("Step 3: The Sufficient Conditions (Contractive Inequalities)")
    print("A unique FGH-tripled fixed point exists if there are non-negative real numbers")
    print(f"({a1}, {a2}, {a3}), ({b1}, {b2}, {b3}), and ({g1}, {g2}, {g3}) such that for any points")
    print("(x₁, x₂ in X), (y₁, y₂, v₁, v₂ in Y), and (z₁, z₂ in Z), the following inequalities hold:\n")
    
    # Printing each inequality with symbolic numbers
    print(f"1. For function F:")
    print(f"   d_X(F(x₁, y₁, z₁), F(x₂, y₂, z₂))  <=  {a1}*d_X(x₁, x₂) + {a2}*d_Y(y₁, y₂) + {a3}*d_Z(z₁, z₂)\n")
    
    print(f"2. For function G:")
    print(f"   d_Y(G(y₁, x₁, v₁), G(y₂, x₂, v₂))  <=  {b1}*d_Y(y₁, y₂) + {b2}*d_X(x₁, x₂) + {b3}*d_Y(v₁, v₂)\n")

    print(f"3. For function H:")
    print(f"   d_Z(H(z₁, y₁, x₁), H(z₂, y₂, x₂))  <=  {g1}*d_Z(z₁, z₂) + {g2}*d_Y(y₁, y₂) + {g3}*d_X(x₁, x₂)\n")
    
    # Step 4: State the final condition for the constants
    print("Step 4: The Final Condition on the Constants")
    print("The collection of constants above must ensure that the overall mapping on the")
    print("product space X * Y * Z is a contraction. This is satisfied if the following condition holds.\n")
    print("First, we define three aggregate coefficients based on the constants from Step 3:")
    print(f"   C_x = {a1} + {b2} + {g3}")
    print(f"   C_y = {a2} + {b1} + {b3} + {g2}")
    print(f"   C_z = {a3} + {g1}\n")
    
    print("The final condition is that the largest of these coefficients must be strictly less than 1.")
    print("This is the crucial final equation that guarantees a unique solution:\n")
    
    final_equation_lhs = f"k = max({a1} + {b2} + {g3}, {a2} + {b1} + {b3} + {g2}, {a3} + {g1})"
    
    print(f"   {final_equation_lhs} < 1")
    

if __name__ == '__main__':
    explain_fgh_tripled_fixed_point_conditions()
