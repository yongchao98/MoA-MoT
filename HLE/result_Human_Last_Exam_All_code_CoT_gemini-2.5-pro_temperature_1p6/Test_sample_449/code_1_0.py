import math

# Plan:
# 1. Define the starting point and calculate its norm.
x0 = (3000, 4000)
norm_x0 = math.sqrt(x0[0]**2 + x0[1]**2)

# 2. Define the constants needed for the asymptotic formula of the potential kernel a(x).
pi = math.pi
gamma = math.e TBD # math.gamma is the gamma function, not the constant. Use the numerical value.
gamma_const = 0.5772156649
ln2 = math.log(2)

# 3. Define the value of the potential kernel at (1,0). This is our reference point.
a_1_0 = 1.0

# 4. Calculate the potential kernel at the starting point x0 using the asymptotic formula.
# a(x) approx= (2/pi) * (ln(|x|) + gamma + 1.5*ln(2))
ln_norm_x0 = math.log(norm_x0)
a_x0 = (2 / pi) * (ln_norm_x0 + gamma_const + 1.5 * ln2)

# 5. Calculate the probability of never hitting the neighbors of the origin.
# P(never hit) = 1 - a(1,0) / a(x0)
prob_never_hit = 1 - a_1_0 / a_x0

# 6. Print the results step-by-step as per the thinking process.
print("The problem is to find the probability that a specific random walk starting at (3000, 4000) never hits the four neighbours of the origin.")
print("This probability P is given by the formula:")
print("P = 1 - a(1,0) / a(3000, 4000)\n")
print("Where a(x) is the potential kernel of the 2D simple random walk.")
print("We use the exact value a(1,0) = 1 and an asymptotic approximation for a(3000, 4000).\n")
print("--- Calculation ---")
print(f"Starting point x0 = {x0}, norm |x0| = {norm_x0}")
print(f"Value of potential kernel at (1,0): a(1,0) = {a_1_0}\n")
print("Asymptotic formula for a(x): a(x) ~= (2/pi) * (ln|x| + gamma + 1.5*ln(2))")
print(f"Constants: pi = {pi:.5f}, gamma = {gamma_const:.5f}, ln(2) = {ln2:.5f}")
print(f"ln|x0| = ln({norm_x0}) = {ln_norm_x0:.5f}\n")
print(f"a(3000, 4000) ~= (2/{pi:.5f}) * ({ln_norm_x0:.5f} + {gamma_const:.5f} + 1.5*{ln2:.5f})")
print(f"a(3000, 4000) ~= {a_x0:.5f}\n")
print("Final Probability Calculation:")
print(f"P = 1 - {a_1_0} / {a_x0:.5f}")
print(f"P = 1 - {a_1_0/a_x0:.5f}")
print(f"P ~= {prob_never_hit:.5f}\n")
print(f"The approximate probability with two significant digits is {prob_never_hit:.2f}")
