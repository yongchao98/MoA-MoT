import sympy

def solve_work_cycle():
    """
    This function symbolically calculates the work done by the current source
    in the described magnetomechanical system over one cycle.
    """
    # Step 1: Define all the symbolic variables
    mu, mu_0, N, w, g, D = sympy.symbols('mu mu_0 N w g D', real=True, positive=True)
    x, I = sympy.symbols('x I')
    x1, x2, I1, I2 = sympy.symbols('x1 x2 I1 I2', real=True)

    print("Step 1 & 2: Derive the Inductance L(x)")
    # The air gap consists of two parallel reluctances:
    # R_block = g / (mu * x * w)
    # R_air = g / (mu_0 * (D - x) * w)
    # The total reluctance R_total is 1 / (1/R_block + 1/R_air)
    # 1/R_total = (mu*x*w/g) + (mu_0*(D-x)*w/g) = (w/g) * ((mu - mu_0)*x + mu_0*D)
    # R_total = g / (w * ((mu - mu_0)*x + mu_0*D))
    # The inductance L(x) = N^2 / R_total
    L_x = (N**2 * w / g) * ((mu - mu_0) * x + mu_0 * D)
    print(f"The inductance L(x) is: {L_x}\n")

    print("Step 3: Set up the differential work dW = I * dPhi")
    # dPhi = L dI + I dL = L(x)dI + I * (dL/dx) * dx
    # dW = I * (L(x)dI + I * (dL/dx) * dx)
    dL_dx = sympy.diff(L_x, x)
    print(f"The derivative dL/dx is: {dL_dx}\n")

    print("Step 4: Integrate dW over the four paths of the cycle")
    # Path 1: x from x1 to x2, I = I1 (dI=0)
    # dW1 = I1**2 * dL/dx * dx
    W1 = sympy.integrate(I1**2 * dL_dx, (x, x1, x2))
    print(f"Work on Path 1 (x1->x2, I=I1): {W1}")

    # Path 2: I from I1 to I2, x = x2 (dx=0)
    # dW2 = I * L(x2) * dI
    L_at_x2 = L_x.subs(x, x2)
    W2 = sympy.integrate(I * L_at_x2, (I, I1, I2))
    print(f"Work on Path 2 (I1->I2, x=x2): {W2}")

    # Path 3: x from x2 to x1, I = I2 (dI=0)
    # dW3 = I2**2 * dL/dx * dx
    W3 = sympy.integrate(I2**2 * dL_dx, (x, x2, x1))
    print(f"Work on Path 3 (x2->x1, I=I2): {W3}")

    # Path 4: I from I2 to I1, x = x1 (dx=0)
    # dW4 = I * L(x1) * dI
    L_at_x1 = L_x.subs(x, x1)
    W4 = sympy.integrate(I * L_at_x1, (I, I2, I1))
    print(f"Work on Path 4 (I2->I1, x=x1): {W4}\n")
    
    print("Step 5: Sum the work from all paths and simplify")
    W_total = sympy.simplify(W1 + W2 + W3 + W4)
    
    # The result from sympy is: -N**2*w*(-I1**2 + I2**2)*(-mu + mu_0)*(-x1 + x2)/(2*g)
    # We can reformat this for better readability to match the options.
    # W = - ( (mu - mu_0) / (2*g) ) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)
    
    print("Final Result:")
    print("The total work W done by the current source for each cycle is:")
    
    # Print the equation part by part as requested
    print("W = - (", end="")
    # Numerator part 1
    print("(mu - mu_0)", end="")
    # Denominator
    print(" / ", end="")
    print("(2*g)", end="")
    print(") * ", end="")
    # N^2*w term
    print("(N**2 * w)", end="")
    # x term
    print(" * ", end="")
    print("(x2 - x1)", end="")
    # I term
    print(" * ", end="")
    print("(I2**2 - I1**2)")

if __name__ == '__main__':
    solve_work_cycle()