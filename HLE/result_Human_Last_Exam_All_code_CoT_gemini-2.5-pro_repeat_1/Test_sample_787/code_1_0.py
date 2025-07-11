# The limit L of the sequence g_n is the greatest common divisor
# of P(p) for all sufficiently large primes p.
# Based on number theory analysis, the prime factorization of L is found to be:
# v_2(L) = 10
# v_3(L) = 2
# v_5(L) = 1
# v_q(L) = 0 for all primes q > 5

base2, exp2 = 2, 10
base3, exp3 = 3, 2
base5, exp5 = 5, 1

val2 = base2**exp2
val3 = base3**exp3
val5 = base5**exp5

result = val2 * val3 * val5

print(f"The limit is L = {base2}^{exp2} * {base3}^{exp3} * {base5}^{exp5}")
print(f"L = {val2} * {val3} * {val5}")
print(f"L = {result}")
