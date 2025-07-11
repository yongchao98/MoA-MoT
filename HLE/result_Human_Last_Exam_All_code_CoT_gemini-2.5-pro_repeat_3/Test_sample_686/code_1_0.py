def derive_susceptibility():
    """
    This function prints the step-by-step derivation of the magnetic susceptibility
    for an Ising model on a sparse random graph.
    """
    print("Derivation of the magnetic susceptibility chi:")
    print("-" * 50)

    # Step 1: Start with the given formula for chi
    print("Step 1: The susceptibility is given by the formula:")
    print("chi = beta * sum_{l=1 to inf} [c * (c-1)^(l-1) * C_l]\n")

    # Step 2: Substitute the expression for C_l
    print("Step 2: We use the relation C_l = (1/beta) * d<sigma_0>/dB_l.")
    print("Substituting this into the formula for chi gives:")
    print("chi = beta * sum_{l=1 to inf} [c * (c-1)^(l-1) * (1/beta) * d<sigma_0>/dB_l]")
    print("chi = sum_{l=1 to inf} [c * (c-1)^(l-1) * d<sigma_0>/dB_l]\n")

    # Step 3: Evaluate the derivative d<sigma_0>/dB_l
    print("Step 3: The derivative d<sigma_0>/dB_l represents the response of the magnetization at site 0")
    print("to a field change at site l. On a tree-like graph, this perturbation propagates along the")
    print("unique path of length l, being multiplied by a transmission factor T at each step.")
    print("The derivative evaluates to:")
    print("d<sigma_0>/dB_l = beta * (1 - m_0^2) * T^l")
    print("where m_0 is the magnetization <sigma_0> and T is the transmission factor.\n")

    # Step 4: Plug the derivative back into the sum
    print("Step 4: Substituting the expression for the derivative into the sum for chi:")
    print("chi = sum_{l=1 to inf} [c * (c-1)^(l-1) * beta * (1 - m_0^2) * T^l]")
    print("We can factor out the terms that do not depend on l:")
    print("chi = beta * c * (1 - m_0^2) * sum_{l=1 to inf} [(c-1)^(l-1) * T^l]\n")

    # Step 5: Sum the geometric series
    print("Step 5: The sum can be identified as a geometric series:")
    print("Sum = T * sum_{l=1 to inf} [(c-1)*T]^(l-1)  (let k = l-1)")
    print("    = T * sum_{k=0 to inf} [(c-1)*T]^k")
    print("This series converges to 1 / (1 - r) with r = (c-1)*T, for |r| < 1.")
    print("So, the sum evaluates to: T / (1 - (c-1)*T)\n")

    # Step 6: Substitute the sum's result
    print("Step 6: Plugging the result of the summation back into the expression for chi:")
    print("chi = beta * c * (1 - m_0^2) * T / (1 - (c-1)*T)\n")

    # Step 7: Use the given constant N for simplification
    print("Step 7: We are given the constant N = beta * c * (1 - m_0^2) / (c-1).")
    print("Rearranging this gives: beta * c * (1 - m_0^2) = N * (c-1).")
    print("Substituting this into our expression for chi:")
    print("chi = N * (c-1) * T / (1 - (c-1)*T)\n")
    
    # Step 8: Final algebraic manipulation
    print("Step 8: For a more compact expression, we perform an algebraic manipulation:")
    print("chi = N * [ ((c-1)*T) / (1 - (c-1)*T) ]")
    print("    = N * [ (1 - (1 - (c-1)*T)) / (1 - (c-1)*T) - 1]")
    print("    = N * [ 1 / (1 - (c-1)*T) - 1 ]\n")

    # Final Answer
    print("-" * 50)
    print("The final expression for the susceptibility is:")
    print("chi = N * (1 / (1 - (c - 1) * T) - 1)")

if __name__ == '__main__':
    derive_susceptibility()