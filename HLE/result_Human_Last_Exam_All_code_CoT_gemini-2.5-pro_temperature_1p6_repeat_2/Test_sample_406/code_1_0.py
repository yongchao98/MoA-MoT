import math

def find_fgh_tripled_fixed_point():
    """
    This function demonstrates the existence of an FGH-tripled fixed point
    by defining a system of functions that satisfies the contraction mapping principle
    and finding the fixed point using iteration.
    """
    
    # 1. DEFINE A SYSTEM SATISFYING THE CONDITIONS
    # We choose simple linear functions for F, G, H over the real numbers (a complete metric space).
    # F(x, y, z) = a1*x + a2*y + a3*z + c1
    # G(y, x, y) = b1*y + b2*x + c2  (Note: the argument 'y' is used twice)
    # H(z, y, x) = g1*z + g2*y + g3*x + c3
    # The absolute values of a_i, b_i, g_i are the Lipschitz constants (α_i, β_i, γ_i).

    # Coefficients for the functions
    a1, a2, a3 = 0.1, 0.1, 0.2
    b1, b2 = 0.2, 0.1
    g1, g2, g3 = 0.1, 0.2, 0.1
    
    # Constant terms
    c1, c2, c3 = 5.0, 2.0, 7.0
    
    # The functions
    F = lambda x, y, z: a1*x + a2*y + a3*z + c1
    G = lambda y, x, _: b1*y + b2*x + c2 # Note: signature adapted to G(y, x, y)
    H = lambda z, y, x: g1*z + g2*y + g3*x + c3

    print("--- Conditions for FGH-Tripled Fixed Point ---")
    print("We have defined a system of functions F, G, H with the following coefficients:")
    print(f"F(x,y,z) = {a1}x + {a2}y + {a3}z + {c1}")
    print(f"G(y,x,y) = {b2}x + {b1}y + {c2}")
    print(f"H(z,y,x) = {g3}x + {g2}y + {g1}z + {c3}\n")
    
    # 2. VERIFY THE CONTRACTION CONDITION
    # A sufficient condition for a unique fixed point to exist is:
    # max( (α1+β2+γ3), (α2+β1+γ2), (α3+γ1) ) < 1
    # where α, β, γ are the absolute values of the coefficients.
    
    k_x = abs(a1) + abs(b2) + abs(g3)
    k_y = abs(a2) + abs(b1) + abs(g2)
    k_z = abs(a3) + abs(g1)
    
    k = max(k_x, k_y, k_z)
    
    print("--- Verifying the Contraction Condition ---")
    print(f"Sum of coefficients for x-dependency (α1+β2+γ3): {k_x:.2f}")
    print(f"Sum of coefficients for y-dependency (α2+β1+γ2): {k_y:.2f}")
    print(f"Sum of coefficients for z-dependency (α3+γ1):    {k_z:.2f}")
    print(f"The contraction constant k = max({k_x:.2f}, {k_y:.2f}, {k_z:.2f}) = {k:.2f}")
    
    if k < 1:
        print(f"Condition met: k = {k:.2f} < 1. A unique fixed point exists.\n")
    else:
        print(f"Condition NOT met: k = {k:.2f} >= 1. Existence is not guaranteed by this theorem.\n")
        return

    # 3. FIND THE FIXED POINT BY ITERATION
    # The Banach theorem guarantees that starting from any point and iterating
    # (x,y,z) -> (F(x,y,z), G(y,x,y), H(z,y,x)) will converge to the fixed point.
    
    x, y, z = 0.0, 0.0, 0.0 # Initial guess
    iterations = 50
    
    print(f"--- Finding the Fixed Point via Iteration (starting from (0,0,0)) ---")
    for i in range(iterations):
        x_next = F(x, y, z)
        y_next = G(y, x, y)
        z_next = H(z, y, x)
        x, y, z = x_next, y_next, z_next
        if i < 5 or (i + 1) % 10 == 0:
            print(f"Iteration {i+1:2d}: (x, y, z) = ({x:.4f}, {y:.4f}, {z:.4f})")
    
    print("\n--- The FGH-Tripled Fixed Point ---")
    print(f"The iterated sequence converges to the fixed point (x, y, z):")
    print(f"x = {x:.8f}")
    print(f"y = {y:.8f}")
    print(f"z = {z:.8f}\n")
    
    print("--- Verification of the Final Equation ---")
    print("This means the following equations hold:")
    print(f"{x:.4f} = {a1}*({x:.4f}) + {a2}*({y:.4f}) + {a3}*({z:.4f}) + {c1}  (which equals {F(x,y,z):.4f})")
    print(f"{y:.4f} = {b2}*({x:.4f}) + {b1}*({y:.4f}) + {c2}  (which equals {G(y,x,y):.4f})")
    print(f"{z:.4f} = {g3}*({x:.4f}) + {g2}*({y:.4f}) + {g1}*({z:.4f}) + {c3}  (which equals {H(z,y,x):.4f})")

if __name__ == '__main__':
    find_fgh_tripled_fixed_point()