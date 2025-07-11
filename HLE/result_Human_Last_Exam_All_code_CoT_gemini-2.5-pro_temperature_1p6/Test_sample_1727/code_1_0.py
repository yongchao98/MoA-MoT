import numpy as np

def calculate_chi():
    """
    This function calculates the stable amplitude factor chi for the given tectonic plate problem.
    
    The analytical solution for the amplitude factor chi is:
    chi = (-kh * cosh(kh)) / (cosh(kh)^2 + kh^2)
    
    Since no specific values for k and h are provided, this function calculates chi
    for the representative case where the dimensionless parameter kh = 1.
    """
    
    # Define the dimensionless parameter kh
    kh = 1.0
    
    # Calculate the terms in the equation
    cosh_kh = np.cosh(kh)
    cosh_sq_kh = cosh_kh**2
    kh_sq = kh**2
    
    numerator = -kh * cosh_kh
    denominator = cosh_sq_kh + kh_sq
    
    # Calculate chi
    chi = numerator / denominator
    
    # Print the results step-by-step
    print("The stable amplitude factor is given by the equation:")
    print("chi = [-kh * cosh(kh)] / [cosh(kh)^2 + kh^2]\n")
    
    print(f"For the representative case where kh = {kh}:\n")
    
    print("--- Numerator Calculation ---")
    print(f"The value of -kh is: {-kh}")
    print(f"The value of cosh(kh) is: {cosh_kh:.6f}")
    print(f"Numerator = {-kh} * {cosh_kh:.6f} = {numerator:.6f}\n")
    
    print("--- Denominator Calculation ---")
    print(f"The value of cosh(kh)^2 is: {cosh_sq_kh:.6f}")
    print(f"The value of kh^2 is: {kh_sq:.6f}")
    print(f"Denominator = {cosh_sq_kh:.6f} + {kh_sq:.6f} = {denominator:.6f}\n")
    
    print("--- Final Calculation ---")
    print(f"chi = {numerator:.6f} / {denominator:.6f}")
    print(f"The final value for the stable amplitude factor chi is: {chi:.6f}")

if __name__ == "__main__":
    calculate_chi()