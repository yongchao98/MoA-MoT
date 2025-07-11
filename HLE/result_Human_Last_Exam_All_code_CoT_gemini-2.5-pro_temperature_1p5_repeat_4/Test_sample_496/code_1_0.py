def F_val(N):
    """Computes sum of coefficients of 1/(1-t^4)^2 up to degree N."""
    if N < 0:
        return 0
    num_terms = N // 4 + 1
    # Sum of arithmetic series m+1 for m from 0 to num_terms-1
    return num_terms * (1 + num_terms) // 2

# Calculate the required F values
F100 = F_val(100)
F98 = F_val(98)
F97 = F_val(97)
F95 = F_val(95)
F91 = F_val(91)

# The total rank equation is:
# Total Rank = F(100) + F(98) + 2*F(97) + F(95) - F(91)
total_rank = F100 + F98 + 2 * F97 + F95 - F91

print("The calculation is based on the formula:")
print("Total Rank = F(100) + F(98) + 2*F(97) + F(95) - F(91)")
print("where F(N) is the sum of ranks of H*(BSO(4) x BSO(4)) up to degree N.")
print("\nIndividual components:")
print(f"F(100) = {F100}")
print(f"F(98) = {F98}")
print(f"F(97) = {F97}")
print(f"F(95) = {F95}")
print(f"F(91) = {F91}")
print("\nFinal equation:")
print(f"Total Rank = {F100} + {F98} + 2*{F97} + {F95} - {F91}")

print("\nResult:")
print(f"The total rank is {total_rank}")
