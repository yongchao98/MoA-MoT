import numpy as np

def solve_lindhard_limit():
    """
    Calculates and explains the numerical value of the dimensionless Lindhard function
    in the static, long-wavelength limit.
    """
    print("This script determines the numerical value of the Lindhard polarization function")
    print("in the static (ω=0) and long-wavelength (q->0) limit for a 3D electron gas at T=0.")
    print("-" * 60)

    print("Step 1: The physical polarization function Π(q,ω).")
    print("The full function, Π(q->0, ω=0), evaluates to -g(ε_F), the negative density of states")
    print("at the Fermi level. This value is not a universal constant.")
    print("\nWe therefore solve for the value of the DIMENSIONLESS Lindhard function, F(x),")
    print("which is defined by the relation: Π(q, 0) = -g(ε_F) * F(q/2k_F).")
    print("-" * 60)

    print("Step 2: The dimensionless static Lindhard function F(x).")
    print("The function is given by: F(x) = 1/2 + ((1 - x^2) / (4*x)) * ln| (1+x)/(1-x) |")
    print("We need to find its limit as momentum q -> 0, which means x -> 0.")
    print("-" * 60)

    print("Step 3: Analytical calculation of the limit.")
    print("To find the limit of F(x) as x approaches 0, we use the Taylor expansion for the logarithm:")
    print("ln|(1+x)/(1-x)| ≈ 2*x + (2/3)*x^3 + O(x^5)")
    print("\nSubstituting this into the equation for F(x):")
    # The numbers in the following equation are 1/2, 1, 4, 2, 0
    one_half = 0.5
    one = 1.0
    four = 4.0
    two = 2.0
    zero = 0.0
    
    print(f"F(x) ≈ {one_half} + (({one} - x^2) / ({four}*x)) * ({two}*x + ...)")
    print(f"F(x) ≈ {one_half} + ({one} - x^2) / {two}")
    print(f"F(x) ≈ {one_half} + {one_half} - x^2/{two}")
    print(f"F(x) ≈ {one} - x^2/{two}")
    print("\nA more accurate expansion F(x) ≈ 1 - x^2/3 shows the behavior near zero.")
    print("\nIn both cases, as x approaches 0, the function F(x) approaches a final value.")
    final_value = 1
    print(f"The equation for the limit is: Lim(x -> {int(zero)}) F(x) = {final_value}")
    print("-" * 60)
    
    print("Step 4: Numerical confirmation.")
    # We define the function F(x) for numerical evaluation
    def dimensionless_lindhard_static(x):
      if abs(x) < 1e-9:
          # Use Taylor expansion for stability near zero
          return 1.0 - (x**2)/3.0
      # Use np.log1p for numerical stability of log(1+y)
      log_term = np.log1p(x) - np.log1p(-x)
      return 0.5 + (1.0 - x**2) / (4.0 * x) * log_term

    x_small = 1e-7
    numerical_result = dimensionless_lindhard_static(x_small)
    print(f"Evaluating F(x) for a very small x (e.g., x = {x_small}):")
    print(f"F({x_small}) = {numerical_result:.15f}")
    print("The numerical result confirms the analytical limit of 1.")
    print("-" * 60)

# Run the solution
solve_lindhard_limit()