import sys

# Disabling the creation of __pycache__
sys.dont_write_bytecode = True

def solve_conormal_space():
    """
    This function determines the conormal space for R(sigma)f.

    The problem asks for the conormal space of u = R(sigma)f, where R(sigma)
    is the resolvent of the wave operator on a Schwarzschild background, and
    f is a function in the conormal space A^{2+alpha}(X).

    Let s_in be the order of the conormal space for f.
    Let s_out be the order of the conormal space for u.

    1. The function f is in A^{s_in}(X) where s_in = 2 + alpha. This implies
       f has an asymptotic decay of O(r^{-s_in}) for large radius r.

    2. The function u = R(sigma)f solves the equation (Box_g - sigma^2)u = f and
       must satisfy the outgoing radiation condition.

    3. The solution u can be seen as a sum of a particular solution (driven by f)
       and a homogeneous solution (to satisfy boundary conditions).
       A particular solution u_p can be found that decays at the same rate as f,
       i.e., O(r^{-s_in}) = O(r^{-(2+alpha)}).

    4. The outgoing radiation condition requires that the solution u has the asymptotic
       form A(omega) * e^{i*sigma*r} / r for some angular function A(omega). This
       term decays like O(r^{-1}).

    5. The overall asymptotic behavior of u is determined by the slowest-decaying part.
       We compare the decay of the particular part, r^{-(2+alpha)}, with the decay
       of the radiating part, r^{-1}.

    6. Since alpha is in (0, 1), the exponent 2 + alpha is strictly greater than 1.
       Therefore, r^{-1} decays slower than r^{-(2+alpha)}.

    7. The dominant term in the asymptotic expansion of u is proportional to r^{-1}.
       In the compactified coordinate rho = 1/r, this corresponds to rho^1.

    8. Thus, the solution u = R(sigma)f belongs to the conormal space A^1(X).
       The output conormal order is s_out = 1.
    """

    # Input conormal order for f. alpha is symbolic.
    s_in_exponent_val_1 = 2
    s_in_exponent_symbol = "alpha"
    s_in = f"{s_in_exponent_val_1} + {s_in_exponent_symbol}"

    # Output conormal order for u = R(sigma)f, based on the reasoning above.
    s_out = 1

    print(f"Let f be in the conormal space A^s(X) with s = {s_in_exponent_val_1} + {s_in_exponent_symbol}.")
    print(f"The resolvent R(sigma)f is the unique outgoing solution to (Box_g - sigma^2)u = f.")
    print(f"The asymptotic behavior of u is determined by the radiation condition, which imposes a decay of O(r^-1).")
    print(f"This corresponds to a conormal order of 1.")
    print(f"Thus, R(sigma)f belongs to the conormal space A^s(X) with the new order s.")
    print(f"Final Equation: s = {s_out}")


if __name__ == '__main__':
    solve_conormal_space()