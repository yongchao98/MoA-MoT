import math

def solve():
    """
    Calculates the value of -u(0,1)/2 using a second-order Taylor expansion for u(0,t).
    """

    # Step 1 & 2: Evaluate u(0,0) and u_t(0,0) from initial conditions.
    # At x=0: tanh(0) = 0, sech(0) = 1, exp(0) = 1
    # u(x,0) = -2 + (1-tanh(x))/(exp(x)+1)
    u_00 = -2 + (1 - math.tanh(0)) / (math.exp(0) + 1)
    # du/dt(x,0) = 1/4 * (tanh(x)-1) * sech(x/2)^2 * (tanh(x)-sech(x)-2)
    # sech(0) = 1
    sech_0_over_2 = 1 / math.cosh(0)
    u_t_00 = 0.25 * (math.tanh(0) - 1) * (sech_0_over_2**2) * (math.tanh(0) - 1 / math.cosh(0) - 2)

    # Step 3: Calculate spatial derivatives u_x(0,0) and u_xx(0,0).
    # These are derived by differentiating u(x,0) and evaluating at x=0.
    # Let f(x) = 1 - tanh(x), g(x) = (e^x+1)^-1. u(x,0) = -2 + f(x)g(x).
    # At x=0:
    # f(0) = 1, f'(0) = -sech^2(0) = -1, f''(0) = 2sech^2(0)tanh(0)=0
    # g(0)=1/2, g'(0)=-e^0/(e^0+1)^2=-1/4, g''(0)=[-e^x(e^x+1)^2+2e^(2x)(e^x+1)]/(e^x+1)^4=0
    # u_x(0,0) = f'(0)g(0) + f(0)g'(0) = (-1)(1/2) + (1)(-1/4) = -0.75
    u_x_00 = -0.75
    # u_xx(0,0) = f''(0)g(0)+2f'(0)g'(0)+f(0)g''(0)=0+2(-1)(-1/4)+0 = 0.5
    u_xx_00 = 0.5

    # Step 4: Use the PDE to find u_tt(0,0)
    # PDE: u_t + 1/8*u_tt + u*u_x - 1/8*u_xx - (u-1)u(u+2) = 0
    # 1/8*u_tt = -u_t - u*u_x + 1/8*u_xx + (u-1)u(u+2)
    P_u_00 = (u_00 - 1) * u_00 * (u_00 + 2)
    u_tt_00_over_8 = -u_t_00 - u_00 * u_x_00 + (1/8) * u_xx_00 + P_u_00
    u_tt_00 = 8 * u_tt_00_over_8

    # Step 5: Approximate u(0,1) using the Taylor expansion
    # u(0,1) approx u(0,0) + u_t(0,0)*t + u_tt(0,0)/2 * t^2 for t=1
    t = 1
    u_01_approx = u_00 + u_t_00 * t + u_tt_00 / 2 * (t**2)

    # Step 6: Calculate the final result
    result = -u_01_approx / 2

    # Print the step-by-step calculation
    print("Step 1: Calculate u(0,0) from the initial condition for u(x,0)")
    print(f"u(0,0) = -2 + (1 - tanh(0)) / (exp(0) + 1) = {u_00}")
    
    print("\nStep 2: Calculate du/dt(0,0) from the initial condition for du/dt(x,0)")
    print(f"du/dt(0,0) = 1/4 * (tanh(0)-1) * sech(0)^2 * (tanh(0)-sech(0)-2) = {u_t_00}")
    
    print("\nStep 3: Calculate spatial derivatives at (0,0)")
    print(f"du/dx(0,0) = {u_x_00}")
    print(f"d^2u/dx^2(0,0) = {u_xx_00}")
    
    print("\nStep 4: Use the PDE to calculate d^2u/dt^2(0,0)")
    print(f"1/8 * d^2u/dt^2(0,0) = -({u_t_00}) - ({u_00})*({u_x_00}) + 1/8*({u_xx_00}) + ({u_00}-1)*({u_00})*({u_00}+2)")
    print(f"1/8 * d^2u/dt^2(0,0) = {-u_t_00 - u_00 * u_x_00 + (1/8) * u_xx_00 + P_u_00}")
    print(f"d^2u/dt^2(0,0) = {u_tt_00}")
    
    print("\nStep 5: Approximate u(0,1) with the Taylor series")
    print(f"u(0,1) ≈ u(0,0) + du/dt(0,0) * (1) + 1/2 * d^2u/dt^2(0,0) * (1)^2")
    print(f"u(0,1) ≈ {u_00} + {u_t_00} + {u_tt_00}/2")
    print(f"u(0,1) ≈ {u_01_approx}")
    
    print("\nStep 6: Calculate the final quantity")
    print(f"The value of -u(0,1)/2 is -({u_01_approx})/2 = {result}")


solve()