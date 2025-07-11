# The problem reduces to calculating 2 * f(2), where f(x) is the second term in the integrand.
# As derived in the explanation, the value of the integral is 2^(15/16).
# Here we calculate this value.

val1 = 2
val2_exp_num = -1
val2_exp_den = 16
val2 = 2**(val2_exp_num / val2_exp_den)

result_exp_num = 15
result_exp_den = 16
result_val = 2**(result_exp_num / result_exp_den)

# The final equation is Value = 2 * 2^(-1/16) = 2^(15/16).
# We print out each number in this final equation as requested.
print("Final Equation:")
print(f"{val1} * 2^({val2_exp_num}/{val2_exp_den}) = 2^({result_exp_num}/{result_exp_den})")
print("\nNumerical Value:")
print(result_val)