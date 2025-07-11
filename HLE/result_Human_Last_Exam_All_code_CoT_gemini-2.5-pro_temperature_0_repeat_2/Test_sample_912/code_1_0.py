import sympy as sp

def solve_work_done():
    """
    This function symbolically calculates the work done by the current source
    over the specified rectangular cycle in the I-x plane.
    """
    # Define all symbols used in the problem.
    # We assume all are positive real numbers.
    mu_0, mu, N, w, g, D, x = sp.symbols('mu_0 mu N w g D x', real=True, positive=True)
    I, I1, I2, x1, x2 = sp.symbols('I I1 I2 x1 x2', real=True, positive=True)

    # Step 1: Model the Magnetic Circuit and find Reluctance R(x)
    # The air gap is composed of two materials in parallel:
    # 1. The magnetic block: area = x*w, length = g, permeability = mu
    # 2. The air: area = (D-x)*w, length = g, permeability = mu_0
    # The reluctance of a path is R = length / (permeability * area).
    R_block = g / (mu * x * w)
    R_air = g / (mu_0 * (D - x) * w)

    # The total reluctance of parallel paths is R_total = 1 / (1/R1 + 1/R2)
    R_total = 1 / (1/R_block + 1/R_air)

    # Step 2: Determine the Inductance L(x)
    # The inductance is given by L = N^2 / R_total.
    L = sp.simplify(N**2 / R_total)
    # L(x) = (N**2 * w / g) * (mu_0*D + (mu - mu_0)*x)

    # Step 3: Calculate the work done over the cycle.
    # The incremental work done by the source is dW = I * d(Phi), where Phi = L(x)*I.
    # dW = I * d(L*I) = I * (L*dI + I*(dL/dx)*dx)
    # We integrate this over the four paths of the cycle.

    # Path 1: I = I1 (constant, dI=0), x from x1 to x2
    dLdx = sp.diff(L, x)
    integrand1 = (I**2 * dLdx).subs(I, I1)
    W1 = sp.integrate(integrand1, (x, x1, x2))

    # Path 2: x = x2 (constant, dx=0), I from I1 to I2
    integrand2 = (I * L).subs(x, x2)
    W2 = sp.integrate(integrand2, (I, I1, I2))

    # Path 3: I = I2 (constant, dI=0), x from x2 to x1
    integrand3 = (I**2 * dLdx).subs(I, I2)
    W3 = sp.integrate(integrand3, (x, x2, x1))

    # Path 4: x = x1 (constant, dx=0), I from I2 to I1
    integrand4 = (I * L).subs(x, x1)
    W4 = sp.integrate(integrand4, (I, I2, I1))

    # Step 4: Sum the work from all paths and simplify.
    W_total = sp.simplify(W1 + W2 + W3 + W4)

    # The symbolic result from sympy is:
    # -N**2*w*(I1**2 - I2**2)*(mu - mu_0)*(x1 - x2)/(2*g)
    # We re-arrange this to match the format of the answer choices.
    # W = - ( (mu - mu_0) / (2*g) ) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)

    # Step 5: Print the final equation.
    # We construct the string representing the final formula.
    term1 = "(mu - mu_0)"
    term2 = "(2*g)"
    term3 = "N**2"
    term4 = "w"
    term5 = "(x2 - x1)"
    term6 = "(I2**2 - I1**2)"

    print("The work done by the current source for each cycle is:")
    print(f"W = - ({term1} / {term2}) * {term3} * {term4} * {term5} * {term6}")

if __name__ == '__main__':
    solve_work_done()