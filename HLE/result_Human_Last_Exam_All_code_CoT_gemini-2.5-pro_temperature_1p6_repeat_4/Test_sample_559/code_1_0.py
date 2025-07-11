import math

def solve():
    """
    This function prints the equation of the separatrix for the given system of ODEs.
    The derivation is explained in the text above. The separatrix is found to be
    a simple parabolic curve.
    """
    
    # The equation of the separatrix is d = -u^2.
    # We can write this in the form A*d^B + C*u^D = E
    
    coeff_d = 1
    power_d = 1
    
    coeff_u = 1
    power_u = 2
    
    constant = 0
    
    # Rearranging d = -u^2 to d + u^2 = 0
    # which is 1*d^1 + 1*u^2 = 0
    
    print("The equation of the separatrix is d = -u**2.")
    print("This can be written in a standard form as: d + u**2 = 0.")
    print("\nThe numbers in this equation are:")
    print(f"Term d: coefficient = {coeff_d}, power = {power_d}")
    print(f"Term u: coefficient = {coeff_u}, power = {power_u}")
    print(f"Constant on the right side: {constant}")
    
    # Final representation of the equation showing each number
    print(f"\nThe equation with all numbers explicitly shown is: ({coeff_d})*d**{power_d} + ({coeff_u})*u**{power_u} = {constant}")

solve()