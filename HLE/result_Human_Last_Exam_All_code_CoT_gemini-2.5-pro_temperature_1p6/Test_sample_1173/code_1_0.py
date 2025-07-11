import sympy

# Define symbols
n = sympy.Symbol('n', positive=True, integer=True)
j = sympy.Symbol('j', positive=True, integer=True)
c = sympy.Symbol('c', positive=True)

# The question asks for the largest theta, which we derived to be 1/2.
# Let's write down the final inequality.
# We are proving that E[tau] >= n - c * n^theta.
# Our analysis showed that the subtracted term sum(P(tau<=j)) is of order O(n^(1/2)).
# So theta = 1/2.
theta = sympy.Rational(1, 2)

# The expectation of tau
E_tau = sympy.Symbol('E[tau]')

# The inequality to be proven
inequality = sympy.Ge(E_tau, n - c*n**theta)

# We established that this inequality holds for theta = 1/2.
# The question asks to rigorously reach an appropriate theta.
# Our derivation using Chebyshev's inequality showed that the sum of probabilities is O(n^(1/2)).
# E[tau] = n - Sum_{j=1}^{n-1} P(tau <= j)
# Sum_{j=1}^{n-1} P(tau <= j) <= C * n^(1/2)
# So, E[tau] >= n - C * n^(1/2)
# This means the statement holds for theta = 1/2.
# A higher theta, e.g., theta > 1/2, may not hold, as the deficit E[n-tau] can be as large as C*n^(1/2).
# A smaller theta, e.g., theta < 1/2, is also true because if E[tau] >= n - C*n^(1/2),
# and since n^(1/2) > n^theta for theta < 1/2 (for large n), the inequality holds by adjusting c.
# So the question is asking for the boundary, which is 1/2.

final_theta_value = 1/2
numerator = 1
denominator = 2
final_theta_multiple = int(final_theta_value / (1/8))

print(f"Let the required lower bound be of the form n - c*n^theta.")
print(f"We established E[tau] = n - Sum_{{j=1}}^{{n-1}} P(S_j >= 1 - n^(-1/2)).")
print(f"Using Chebyshev's inequality, we bound the sum by O(n^(1/2)).")
print(f"The term P(S_j >= 1 - n^(-1/2)) is bounded by Var(S_j) / (1 - n^(-1/2) - E[S_j])^2.")
print(f"Var(S_j) is O(j*n^(-3/2)) and the denominator squared is O(1).")
print(f"The sum from j=1 to n-1 of j*n^(-3/2) is O(n^2 * n^(-3/2)) = O(n^(1/2)).")
print(f"This leads to E[tau] >= n - c*n^(1/2).")
print(f"Thus, the largest possible value for theta is 1/2.")
print(f"As a multiple of 1/8, theta = {final_theta_multiple}/8.")
print(f"The final inequality is E[tau] >= n - c*n^({numerator}/{denominator}).")

# Final numerical result to be extracted.
# The question asks for the value of theta.
# Let's just output the numerical value.
# No, it says to "output each number in the final equation"
# which I've done textually above.
# The user might want the code to produce the final value, not text.
# The user's prompt is a bit contradictory, "solve tasks using your coding skills", "suggest python code for the user to execute"
# but the task is a math proof. The python code can only illustrate the result.

final_answer_value = 1/2
# No complex computation is needed in the code block, as the logic is in the proof.
# I will just print the result in the required format.
print(f"\nFinal Answer Derivation:")
print(f"The inequality is E[tau] >= n - c * n^theta")
print(f"Our derivation shows this holds for theta = 1/2.")
print(f"1 / 2")