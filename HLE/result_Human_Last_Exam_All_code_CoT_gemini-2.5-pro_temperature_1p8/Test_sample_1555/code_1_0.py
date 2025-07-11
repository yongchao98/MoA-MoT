import math

# Plan:
# 1. State the mass formula for excitations in the CP(N-1) model.
# 2. Identify the masses of the lightest (k=1) and subsequent (k=2) excitations.
# 3. Form the ratio m_2 / m_1.
# 4. Calculate the limit of this ratio as N -> infinity.

print("Step 1: The mass spectrum formula")
print("In the large-N limit, the mass of the k-th solitonic excitation (m_k) in the CP(N-1) model is:")
print("m_k = m_1 * (sin(pi * k / N) / sin(pi / N))")
print("where m_1 is the mass of the lightest excitation (for k=1).\n")

print("Step 2: Identify the relevant masses")
print("The lightest excitation corresponds to k=1, with mass m_1.")
print("The subsequent higher excitation corresponds to k=2, with mass m_2.")
print("Using the formula for k=2:")
print("m_2 = m_1 * (sin(pi * 2 / N) / sin(pi / N))\n")

print("Step 3: Formulate the mass ratio")
print("The desired mass ratio R is m_2 / m_1.")
print("R = [m_1 * (sin(2 * pi / N) / sin(pi / N))] / m_1")
print("R = sin(2 * pi / N) / sin(pi / N)\n")

print("Step 4: Calculate the asymptotic limit as N -> infinity")
print("We evaluate the limit of R as N approaches infinity.")
print("Let x = pi / N. As N -> oo, x -> 0.")
print("The ratio becomes: R = sin(2x) / sin(x)")
print("Using the trigonometric identity sin(2x) = 2 * sin(x) * cos(x):")
print("R = (2 * sin(x) * cos(x)) / sin(x) = 2 * cos(x)\n")

print("Now, we take the limit as x approaches 0:")
print("Asymptotic Ratio = lim_{x->0} (2 * cos(x))\n")

# Define the numbers for the final equation to be printed
coefficient = 2
limit_point = 0
cos_of_limit_point = math.cos(limit_point)
final_result = coefficient * cos_of_limit_point

print("Final calculation:")
# Output each number in the final equation: 2 * cos(0) = 2 * 1 = 2
print(f"The ratio is {coefficient} * cos({limit_point}) = {coefficient} * {int(cos_of_limit_point)} = {int(final_result)}")
