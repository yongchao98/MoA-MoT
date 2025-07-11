import math

def solve_blowup_condition():
    """
    Analyzes the blow-up conditions for the given system of ODEs and prints the result.

    The system of differential equations is:
    x'(t) = -3*x(t)*y(t)
    y'(t) = -y(t)^2 - x(t) + 1

    We are given the initial condition x(0) > 1 and need to find the values of y(0)
    for which the solution (x(t), y(t)) blows up.
    """

    print("### Analysis of Blow-Up Conditions ###")
    print("\nStep 1: Analyze the case where y(0) <= 0.")
    print("If y(0) <= 0, it can be shown that y(t) will remain negative and decrease, while x(t) increases from x(0) > 1.")
    print("The equation for y'(t) becomes strongly negative, leading to a finite-time blow-up where y(t) -> -infinity.")
    print("Conclusion: For all y(0) <= 0, the solution blows up.\n")

    print("Step 2: Analyze the case where y(0) > 0.")
    print("Solutions that do not blow up must approach an equilibrium point. The only relevant equilibrium is the saddle point at (1, 0).")
    print("A solution converges to (1, 0) only if its initial condition lies on the stable manifold (a specific path leading to the saddle point).")
    print("To find this path, we first find a conserved quantity Psi(x, y) which is constant on any trajectory.")
    print("By finding an integrating factor for the path equation dy/dx, the conserved quantity is found to be:")
    print("Psi(x, y) = (3/2)*x^(-2/3)*y^2 - 3*x^(1/3) - (3/2)*x^(-2/3)\n")

    print("Step 3: Find the equation for the stable manifold.")
    print("The stable manifold corresponds to the trajectory passing through (1, 0).")
    print("We evaluate Psi at (1, 0) to find the constant value for this specific path.")
    psi_at_saddle = -3.0 - 1.5
    print(f"Psi(1, 0) = -3 - 3/2 = {psi_at_saddle}")
    print("So, any initial condition (x(0), y(0)) on the stable manifold must satisfy Psi(x(0), y(0)) = -4.5.")
    print("This gives the equation: (3/2)*x(0)^(-2/3)*y(0)^2 - 3*x(0)^(1/3) - (3/2)*x(0)^(-2/3) = -4.5")
    print("Solving for y(0)^2, we get the equation for the initial conditions that do NOT lead to blow-up:")
    
    # Define the numbers in the final equation
    c1, c2, c3, n1, n2 = 2, 1, 3, 2, 3
    
    print(f"y(0)^2 = {c1}*x(0) + {c2} - {c3}*x(0)^({n1}/{n2})\n")
    
    print("Step 4: Determine the final condition for blow-up.")
    print("Analysis of the vector field shows that only the positive branch of this manifold is stable.")
    print("Therefore, for any given x(0) > 1, there is exactly ONE positive value of y(0) for which the solution does not blow up.")
    print("This value is given by the equation:")
    print(f"y(0) = sqrt({c1}*x(0) + {c2} - {c3}*x(0)^({n1}/{n2}))\n")
    print("For any other value of y(0), the solution trajectory will not converge to the equilibrium and will eventually enter the blow-up region (x > 1, y < 0).")
    
    final_conclusion = "The solution blows up for all real values of y(0) except for the unique positive value sqrt(2*x(0) + 1 - 3*x(0)^(2/3))."
    print("### Conclusion ###")
    print(final_conclusion)


solve_blowup_condition()

# The final answer is a set of values, not a single numerical value.
# The set is all real numbers except for one specific value that depends on x(0).
# We represent this condition as the answer.
final_answer_string = "y(0) != sqrt(2*x(0) + 1 - 3*x(0)^(2/3))"
print(f"\n<<<For a given x(0) > 1, the solution blows up if {final_answer_string}>>>")