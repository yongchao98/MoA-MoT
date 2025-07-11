import sys
import math

print("Analyzing each statement about Python's division semantics:")
print("-" * 60)

# A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3
x, y = 10, 3
is_false_in_py3 = (x / y) != (x // y)
print(f"A) In Python 3: {x} / {y} == {x} // {y} -> {x/y} == {x//y} is {not is_false_in_py3}.")
print(f"   The statement is TRUE because in Py2, integer division made 10/3 -> 3, which equals 10//3 -> 3. In Py3, it's false.")
print("-" * 60)

# B) The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions
lhs_b = 1/2 * 2
rhs_b = 1/2 + 1/2
print(f"B) In Python 3: 1/2 * 2 == 1/2 + 1/2 -> {lhs_b} == {rhs_b} is {lhs_b == rhs_b}.")
print(f"   The statement is TRUE. (In Py2, it would be 0 == 0, also true).")
print("-" * 60)

# C) For integers a, b: a / float(b) has always returned same result as float(a) / b
a, b = 10, 3
lhs_c = a / float(b)
rhs_c = float(a) / b
print(f"C) {a} / float({b}) == float({a}) / {b} -> {lhs_c} == {rhs_c} is {lhs_c == rhs_c}.")
print(f"   The statement is TRUE. Promoting either operand to a float results in float division.")
print("-" * 60)

# D) In Python 3, -7 // 3 and int(-7 / 3) always give same result
lhs_d = -7 // 3
rhs_d = int(-7 / 3)
print(f"D) -7 // 3 == int(-7 / 3) -> {lhs_d} == {rhs_d} is {lhs_d == rhs_d}.")
print(f"   The statement is FALSE. Floor division (-7//3) rounds to -3, while int() truncates (-2.33... to -2).")
print("-" * 60)

# E) For any non-zero integers x, y: x / y == float(x) / float(y) in Python 3
x, y = 10, 3
lhs_e = x / y
rhs_e = float(x) / float(y)
print(f"E) In Python 3: {x} / {y} == float({x}) / float({y}) -> {lhs_e} == {rhs_e} is {lhs_e == rhs_e}.")
print(f"   The statement is TRUE. The '/' operator performs float division by definition in Python 3.")
print("-" * 60)

# F) The expression 5.0 // 2 has always returned a float across all versions
res_f = 5.0 // 2
is_float_f = isinstance(res_f, float)
print(f"F) The result of 5.0 // 2 is {res_f}. Its type is float: {is_float_f}.")
print(f"   The statement is TRUE. If any operand in // is float, the result is float.")
print("-" * 60)

# G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
n, m = 13, 4
lhs_g = (n//m)*m + n%m
rhs_g = n
print(f"G) ({n}//{m})*{m} + {n}%{m} == {n} -> {lhs_g} == {rhs_g} is {lhs_g == rhs_g}.")
print(f"   The statement is TRUE. This is the division algorithm identity.")
print("-" * 60)

# H) In Python 3, -3.0 // 2.0 == -(3.0 // 2.0)
lhs_h = -3.0 // 2.0
rhs_h = -(3.0 // 2.0)
print(f"H) -3.0 // 2.0 == -(3.0 // 2.0) -> {lhs_h} == {rhs_h} is {lhs_h == rhs_h}.")
print(f"   The statement is FALSE.")
print("-" * 60)

# I) For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)
a, b = 6, 10
lhs_i = a / b
rhs_i = a * (1.0/b)
print(f"I) {a} / {b} == {a} * (1.0/{b}) -> {lhs_i} == {rhs_i} is {lhs_i == rhs_i}.")
print(f"   The statement is FALSE due to floating point representation differences.")
print("-" * 60)

# J) The expression sys.maxsize / 2 returns float in Python 3 even when result fits in int
res_j = sys.maxsize / 2
is_float_j = isinstance(res_j, float)
print(f"J) The type of sys.maxsize / 2 is float: {is_float_j}.")
print(f"   The statement is TRUE. The '/' operator always returns a float in Python 3.")
print("-" * 60)

# K) Floor division of negative numbers rounds toward negative infinity in all versions
num, den = -7, 3
lhs_k = num // den
rhs_k = math.floor(num / den)
print(f"K) {num} // {den} == math.floor({num}/{den}) -> {lhs_k} == {rhs_k} is {lhs_k == rhs_k}.")
print(f"   The statement is TRUE. This is the definition of Python's floor division.")
print("-" * 60)

# L) For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ
n = 5
res_l = (n/1) - n
print(f"L) For n={n}, (n/1)-n results in {res_l} (type {type(res_l)}), which equals 0: {res_l == 0}.")
print(f"   The statement is TRUE. The value is numerically 0, but the type is float.")
print("-" * 60)

# M) Given integers a,b,c: if (a/b)/c == a/(b*c) is true numerically, it's also true for their types
a,b,c = 8, 2, 2
type_lhs_m = type((a/b)/c)
type_rhs_m = type(a/(b*c))
print(f"M) The type of ({a}/{b})/{c} is {type_lhs_m} and the type of {a}/({b}*{c}) is {type_rhs_m}.")
print(f"   The types are equal: {type_lhs_m == type_rhs_m}. The statement is TRUE.")
print("-" * 60)

# N) For positive integers x,y: sum([x/y]*y) equals x in Python 3 only if x is divisible by y
x, y = 10, 4 # Counterexample: 10 is not divisible by 4
sum_res_n = sum([x/y] * y)
is_eq_n = (sum_res_n == x)
is_div_n = (x % y == 0)
print(f"N) For x={x}, y={4}: 'sum([{x}/{y}]*{y}) == {x}' is {is_eq_n}.")
print(f"   '...only if x is divisible by y' ('{x}%{y}==0' is {is_div_n}).")
print(f"   The statement is FALSE, as we found a case where the sum is equal to x, but x is not divisible by y.")
print("-" * 60)

# O) The expression (1/3) + (1/3) + (1/3) == 1.0 evaluates to True even with floating-point arithmetic
lhs_o = (1/3) + (1/3) + (1/3)
rhs_o = 1.0
print(f"O) (1/3)+(1/3)+(1/3) == 1.0 -> {lhs_o} == {rhs_o} is {lhs_o == rhs_o}.")
print(f"   The statement is FALSE due to floating-point imprecision.")
print("-" * 60)