import math

def solve_susceptibility():
    """
    This function prints the step-by-step derivation of the magnetic susceptibility
    for an Ising model on a sparse random graph and outputs the final formula.
    """
    
    print("Derivation of the Magnetic Susceptibility (chi)\n")
    
    print("Step 1: Start with the given expression for chi.")
    print("chi = beta * Sum_{l=1 to inf} [c * (c-1)^(l-1) * C_l]")
    print("Using C_l = (1/beta) * d<sigma_0>/dB_l, we get:")
    print("chi = Sum_{l=1 to inf} [c * (c-1)^(l-1) * d<sigma_0>/dB_l]\n")

    print("Step 2: Find the derivative d<sigma_0>/dB_l.")
    print("By linearizing the message passing equations, one can show that the influence of a field at site l on site 0, which is d<sigma_0>/dB_l, is given by:")
    print("d<sigma_0>/dB_l = beta * (1 - m_0^2) * A^l")
    print("where A is the propagation factor for the perturbation.\n")

    print("Step 3: Apply the small field limit (B -> 0).")
    print("In this limit, the magnetization m_0 = 0, and the propagation factor simplifies to A = tanh(beta*J).\n")
    
    print("Step 4: Substitute the simplified derivative back into the sum.")
    print("d<sigma_0>/dB_l at B=0 becomes: beta * (tanh(beta*J))^l")
    print("So, chi = Sum_{l=1 to inf} [c * (c-1)^(l-1) * beta * (tanh(beta*J))^l]")
    print("Factoring out constants: chi = c * beta * tanh(beta*J) * Sum_{l=1 to inf} [(c-1)*tanh(beta*J)]^(l-1)\n")

    print("Step 5: Evaluate the geometric series.")
    print("The sum is a geometric series Sum_{k=0 to inf} x^k = 1/(1-x), with x = (c-1)*tanh(beta*J).")
    print("This converges for x < 1, which is the condition for the paramagnetic phase.\n")

    print("Step 6: State the final expression.")
    print("The final result for the magnetic susceptibility is:\n")
    
    # Print the final equation, highlighting the numbers as requested.
    one = 1
    print(f"chi = (c * beta * tanh(beta * J)) / ({one} - (c - {one}) * tanh(beta * J))")

solve_susceptibility()