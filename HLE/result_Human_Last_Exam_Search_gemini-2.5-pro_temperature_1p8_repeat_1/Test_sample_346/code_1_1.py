import os
os.system('pip install dovekie')

from dovekie import S, K, I
import math

# The expression, transcribed carefully.
# S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))

f_arg = S(S(K(S))(K))(I)
s_i_i = S(I)(I)
g_arg = S(S(K(S))(K))(s_i_i)

part1 = S(I)
part2 = S(I)(S(I)(K(f_arg)))
part3 = g_arg

expr = part1(part2)(part3)

# Function to convert Church numeral to integer
church_to_int = lambda n: n(lambda x: x + 1)(0)

try:
    n = church_to_int(expr)
    if n > 0:
        log_2_n = math.log2(n)
        print(f"The integer n is: {n}")
        print(f"The value of log_2(n) is: {log_2_n}")
    else:
        print(f"The expression evaluated to n = {n}, which is not a positive integer for which log_2 can be calculated.")
except Exception as e:
    print(f"An error occurred during evaluation: {e}")
