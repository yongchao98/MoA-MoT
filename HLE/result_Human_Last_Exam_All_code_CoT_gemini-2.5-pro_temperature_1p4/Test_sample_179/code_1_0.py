import math

def solve_bvp():
    """
    Solves for X_0(T) for the given boundary-value problem.
    """
    # Given parameters
    A = 1e10
    B = 1e-5 - 1
    T = 1e20
    alpha1 = 0

    print("Step 1: Set up the differential equation for X_0(t).")
    print("The equation is X'_0(t) = -(B + 1)X_0(t) + A.")
    
    # The coefficient k = B + 1
    k = B + 1
    
    print(f"With A = {A} and B = {B}, we have k = B + 1 = {k}.")
    print(f"The equation becomes: X'_0(t) + {k}*X_0(t) = {A}\n")
    
    print("Step 2: State the general solution.")
    print("The general solution for this first-order linear ODE is:")
    print("X_0(t) = C * exp(-k*t) + A/k, where C is the integration constant.")
    
    # Calculate the particular/steady-state solution
    X_ss = A / k
    print(f"The particular solution component is A/k = {A} / {k} = {X_ss:.1e}.\n")

    print("Step 3: Apply the boundary condition X_0(0) = X_0(T).")
    print(f"X_0(0) = C * exp(0) + {X_ss:.1e} = C + {X_ss:.1e}")
    print(f"X_0(T) = C * exp(-k*T) + {X_ss:.1e}")
    print("Setting X_0(0) equal to X_0(T) gives:")
    print(f"C + {X_ss:.1e} = C * exp(-k*T) + {X_ss:.1e}")
    print("C = C * exp(-k*T)  =>  C * (1 - exp(-k*T)) = 0\n")

    print("Step 4: Solve for the constant C.")
    # Value of the exponent
    kT = k * T
    print(f"The term k*T = {k} * {T} = {kT:.1e}.")
    # Because kT is a large positive number, exp(-kT) is extremely close to 0.
    # Therefore, (1 - exp(-kT)) is not zero.
    print(f"Since k*T is a large positive number, exp(-k*T) is practically zero.")
    print("The only way for C * (1 - exp(-k*T)) = 0 to be true is if C = 0.")
    print("This means the solution is unique and constant: X_0(t) = A/k.\n")
    
    print("Step 5: Final calculation of X_0(10^20).")
    final_value = X_ss
    print("The value of the solution at any time t, including T, is A/(B+1).")
    print(f"X_0(T) = A / (B + 1)")
    print("Substituting the values:")
    print(f"X_0({T:.0e}) = {A:.0e} / ({B} + 1) = {A:.0e} / {k} = {final_value:.1e}")

solve_bvp()