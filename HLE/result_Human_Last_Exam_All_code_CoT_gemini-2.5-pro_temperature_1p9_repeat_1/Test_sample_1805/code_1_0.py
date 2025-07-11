# The user can change this value to calculate the q-difference quotient for a different integer power.
n = 4

print(f"Calculating nabla_q(T^n) for n = {n}\n")

if not isinstance(n, int) or n < 0:
    print("Error: This script is designed for non-negative integers n.")
# Case n = 0
elif n == 0:
    # nabla_q(T^0) = nabla_q(1) = (1 - 1) / ((q-1)T) = 0
    print("nabla_q(T^0) = 0")

# Case n = 1
elif n == 1:
    # nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
    print("nabla_q(T^1) = 1")

# Case n > 1
else:
    # Build the T^(n-1) part of the string
    power_of_T = n - 1
    if power_of_T == 1:
        t_part = "T"
    else:
        t_part = f"T^{power_of_T}"

    # Build each term of the expanded sum
    terms = []
    for k in range(n):
        # Term for k=0 ([n]_q starts with 1*T^(n-1))
        if k == 0:
            terms.append(t_part)
        # Term for k=1 (second term is q*T^(n-1))
        elif k == 1:
            terms.append(f"q * {t_part}")
        # General term for k>1 (q^k * T^(n-1))
        else:
            terms.append(f"q^{k} * {t_part}")
    
    # Join all terms with " + "
    final_expression = " + ".join(terms)
    
    print(f"nabla_q(T^{n}) = {final_expression}")
