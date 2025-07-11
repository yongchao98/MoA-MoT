import math

def calculate_blowup_time(a):
    """
    Calculates the blow-up time for the linearized equation ∂_t u + (1+t)Δu = 0.

    The analysis shows that for an initial condition whose Fourier transform decays
    like exp(-a^2|k|^2), the solution's L2-norm blows up when t^2 + 2t - a^2 = 0.
    This function computes the positive root of this quadratic equation, which
    corresponds to the finite blow-up time.

    Args:
        a (float): A parameter characterizing the initial condition's spread in
                   Fourier space. A smaller 'a' implies more high-frequency content.
                   Must be positive.
    """
    if a <= 0:
        print("Parameter 'a' must be positive.")
        return

    # The equation for blow-up time T is T^2 + 2*T - a^2 = 0.
    # We solve for the positive root T = (-2 + sqrt(4 + 4*a^2)) / 2 = -1 + sqrt(1 + a^2).
    c = -a**2
    b = 2.0
    
    # Using the quadratic formula for T^2 + b*T + c = 0
    discriminant = b**2 - 4 * 1 * c
    T_blowup = (-b + math.sqrt(discriminant)) / 2.0
    
    print(f"For an initial condition with characteristic parameter a = {a},")
    print("the blow-up is determined by the equation: T^2 + (2.0)*T + ({c}) = 0")
    print(f"The calculated finite blow-up time T is: {T_blowup}")

# Example: Choose a characteristic parameter for the initial data.
# A smaller 'a' means the initial data is more concentrated in high frequencies,
# leading to a faster blow-up. Let's use a=1 as an example.
a_param = 1.0
calculate_blowup_time(a_param)

# Example with more high-frequency content
# a_param_2 = 0.5
# calculate_blowup_time(a_param_2)