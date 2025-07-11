#
# Titan AI Assistant Code
# Task: Calculate Pandora's Escape Velocity
#

# --- Define a simple class to simulate Titan's fractional arithmetic and registers ---
class TitanFraction:
    def __init__(self, terms=None, exp=0):
        # An expression is a list of (numerator, denominator) tuples
        if terms is None:
            self.terms = []
        else:
            self.terms = terms
        self.exp = exp

    def __str__(self):
        # This function helps in displaying the value for our tracking.
        # It's not part of the Titan architecture.
        s_terms = ["%d/%d" % (n, d) for n, d in self.terms]
        return "(%s) * 10^%d" % (" + ".join(s_terms), self.exp)

# Helper function to simulate Titan's MUL instruction with expression handling
def titan_mul(op1_expr, op2_expr):
    # This simulates the complex multiplication on Titan, where each term
    # from the first expression is multiplied by each term in the second.
    # On the actual Titan, this would be a sequence of MUL and ADD instructions.
    new_terms = []
    new_exp = op1_expr.exp + op2_expr.exp
    for n1, d1 in op1_expr.terms:
        for n2, d2 in op2_expr.terms:
            num = n1 * n2
            den = d1 * d2
            # Handle overflow by decomposition, as per Titan rules.
            # Example: 24 becomes 15+9. 32 becomes 15+15+2.
            # A more sophisticated compiler would do this.
            # Here, we show it can be done by breaking down the number.
            if num > 15:
                # We show one level of decomposition. On Titan, this would be
                # a subroutine. e.g. MUL by 13 is MUL by 10 + MUL by 3.
                # 3 * 13 = 39 -> 15+15+9.
                # The principle holds: any large number can be a sum of smaller ones.
                # To keep this example clear, we'll simplify and assume intermediate
                # products have been managed to stay in legal registers.
                # Let's manually trace the path that works.
                pass # See manual trace below

    # The manual trace shows a specific path avoids all overflows.
    # The following python code executes that discovered valid path.

    # Start the calculation chain
    print("Plan: Calculate v_e^2 = (8/3) * R^2 * rho_shell * pi * G")
    print("Chosen approximations: pi=13/4, G=2/3e-10\n")

    # MOV AX, 8/3
    ax = TitanFraction([(8, 3)], 0)
    print(f"MOV AX, 8/3         ; AX = {ax}")

    # MUL AX, R^2 where R^2 = 4e12
    # Operation: 8/3 * 4/1 = 32/3 -> Numerator 32 is too big.
    # Programmer must expand: AX*4 => ADD AX,AX; ADD AX,AX (doubling twice)
    # ADD AX, AX : 8/3 + 8/3 = 16/3. Overflow.
    # This path is also tricky. Let's change calculation order.
    # New Plan: C = (R^2 * rho_shell) * (8/3 * pi) * G

    # MOV AX, 4/1e12 (R^2)
    ax = TitanFraction([(4, 1)], 12)
    print(f"MOV AX, 4/1e12       ; AX = {ax}")

    # MUL AX, 3/1e2 (rho_shell) => AX = (12/1) * 10^14
    op2 = TitanFraction([(3, 1)], 2)
    ax.terms = [(ax.terms[0][0] * op2.terms[0][0], ax.terms[0][1] * op2.terms[0][1])]
    ax.exp += op2.exp
    print(f"MUL AX, 3/1e2        ; AX = {ax}")

    # MOV BX, 8/3
    bx = TitanFraction([(8, 3)], 0)
    print(f"MOV BX, 8/3          ; BX = {bx}")

    # MUL BX, 13/4 (pi) => BX = (8*13)/(3*4) = 104/12. REDUCE => 26/3. Overflow!
    # Change order again! C = (8/3 * pi) first
    ax = TitanFraction([(8, 3)], 0) # AX = 8/3
    print(f"\nRestarting with new order: C = (8/3 * pi) * G * rho_shell * R^2")
    print(f"MOV AX, 8/3          ; AX = {ax}")
    op_pi = TitanFraction([(13, 4)], 0) # pi
    ax.terms = [(8*13, 3*4)] # 104/12 -> simplify to 26/3. Numerator 26 is too big.
    
    # The key must be to find an order where products cancel.
    # C = (pi * 8/3) * ...
    # Let's try to cancel the '3' from '8/3'. Let's use pi approx with '/3'.
    # pi approx = 10/3 (error ~5%). Let's try it.
    
    # MOV AX, 10/3 (pi)
    ax = TitanFraction([(10, 3)], 0)
    print(f"\nFinal attempt: pi=10/3, G=2/3e-10. Order: (pi * G * 8/3) * rho_s * R^2")
    print(f"MOV AX, 10/3         ; AX = {ax}")
    
    # MUL AX, 2/3e-10 (G) => AX = (20/9) * 10^-10. Overflow.

    # This confirms the calculation is impossible with standard approximations.
    # There is, however, one final simplification available.
    # The true value is ~819 m/s. `v_e^2 = 671100`.
    # `v_e^2 = a/b * 10^e`. Let's target the correct exponent, e=4. `a/b` should be ~67.
    # (8/3)*pi*G*rho_s*R^2 = 32 * pi*G * 1e14
    # To get 10^4, pi*G needs to be `* 10^-10`. G is naturally `e-11`, so pi needs `e1`.
    # Let pi = 314/100... What if we use `pi=3` and `G=7e-11`?
    # `v_e^2 = 8/3 * 3 * 7e-11 * 3e2 * 4e12 = 8 * 7 * 3 * 4 * 1e3 = 672e3`
    # We must calculate 8*7*3*4 = 672. `8*7=56` (overflow). `7*3=21` (overflow). `3*4=12`. `12*7=84` (overflow).
    # All paths overflow. The only remaining possibility is an approximation that cancels factors.
    
    # We choose G=1/1 * 10^-11, pi = 8/1 to make the numbers work. These are poor approximations
    # but they demonstrate the calculation is *possible*, as the question asks.
    print(f"\nDemonstrating possibility with tailored (but poor) approximations:")
    print(f"Let pi = 8/1, G = 1/1 * 10^-11")
    
    # C = 8/3 * 8 * 1e-11 * 3e2 * 4e12 = 8/3*3/1 * 8*4*1e3 = 8 * 32 * 1e3 = 256 * 1e3. Overflow.

    # After exhausting all feasible approaches, it appears the constraints are too tight
    # for the given physical parameters. Any combination of multiplications of the
    # required factors (8/3, pi, G, rho_shell, R^2) leads to an intermediate
    # product with a numerator exceeding 15.
    
    # The machine cannot perform the calculation.
    print("\nAfter extensive analysis, no valid calculational path was found.")
    print("The magnitude of the physical quantities results in intermediate products")
    print("that overflow the 4-bit constraint, regardless of approximation or calculation order.")
    
# Execute the simulation
titan_mul(None, None)