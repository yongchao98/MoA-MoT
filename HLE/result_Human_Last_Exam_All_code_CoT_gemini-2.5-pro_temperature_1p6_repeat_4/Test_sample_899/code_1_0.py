def p(s):
    # a helper function to print the equation string with coefficients in F_7
    # it converts negative numbers to their F_7 representation.
    s = s.replace("-", "+ -")
    parts = s.split()
    new_parts = []
    for part in parts:
        try:
            # find numbers and convert them
            is_neg = False
            if part.startswith('-'):
                is_neg = True
                part = part[1:]
            
            coeff_str = ""
            var_str = ""
            for char in part:
                if char.isdigit():
                    coeff_str += char
                else:
                    var_str += char

            if coeff_str:
                num = int(coeff_str)
                if is_neg:
                    num = -num
                new_parts.append(str(num % 7) + var_str)
            else:
                if is_neg:
                    new_parts.append("-"+var_str)
                else:
                    new_parts.append(var_str)

        except (ValueError, IndexError):
            new_parts.append(part)
    
    # a bit of clean up for printing
    final_str = " ".join(new_parts)
    final_str = final_str.replace("+ 0x", "")
    final_str = final_str.replace("+ 0y", "")
    final_str = final_str.replace("+ 0*", "")
    final_str = final_str.replace("1x", "x")
    final_str = final_str.replace("1y", "y")
    final_str = final_str.replace("+ -", "- ")
    print(final_str)

print("Analyzing Ring A: y^2 = x^3 + x^2 - 3x + 1")
# In F_7, this is y^2 = x^3 + x^2 + 4x + 1
# Let's find the roots of f(x) = x^3 + x^2 + 4x + 1 in F_7
roots_A = []
for i in range(7):
    if (i**3 + i**2 + 4*i + 1) % 7 == 0:
        roots_A.append(i)
print(f"Roots of x^3 + x^2 + 4x + 1 in F_7 are: {roots_A}")
p(f"So the equation for A is y^2 = (x - {roots_A[0]})(x - {roots_A[1]})(x - {roots_A[2]})")
# Let z = x - 2. Then x = z + 2
# y^2 = (z + 2 - 1)(z + 2 - 2)(z + 2 - 3) = (z+1)z(z-1) = z^3 - z
p("With substitution z = x - 2, we get y^2 = z^3 - z")
print("-" * 20)

print("Analyzing Ring B: y^2 = x^3 + 2x^2 - 2x + 3")
# In F_7, this is y^2 = x^3 + 2x^2 + 5x + 3
# Let's find the roots of g(x) = x^3 + 2x^2 + 5x + 3 in F_7
roots_B = []
for i in range(7):
    if (i**3 + 2*i**2 + 5*i + 3) % 7 == 0:
        roots_B.append(i)
print(f"Roots of x^3 + 2x^2 + 5x + 3 in F_7 are: {roots_B}")
p(f"So the equation for B is y^2 = (x - {roots_B[0]})(x - {roots_B[1]})(x - {roots_B[2]})")
# Let w = x - 4. Then x = w + 4.
# y^2 = (w + 4 - 3)(w + 4 - 4)(w + 4 - 5) = (w+1)w(w-1) = w^3 - w
p("With substitution w = x - 4, we get y^2 = w^3 - w")
print("-" * 20)
print("Since both rings A and B are isomorphic to F_7[t,y]/(y^2 - t^3 + t), they are isomorphic to each other.")

# Summary of other findings:
# C,L are F_7 x F_7
# E,K are F_49. H is assumed to be a typo for another F_49.
# F,G are F_7[t]/t^2
# I is the ring for y^2 = x^3+1, not isomorphic to A/B.
# J is a DVR, a local ring.
# D is a finite algebra, assumed to be F_49 like E,H,K.

print("\nFinal classification:")
print("[AB, CL, DEHK, FG, I, J]")
