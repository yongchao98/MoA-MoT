import cypari2

# Initialize the PARI/GP environment through cypari2
pari = cypari2.pari_from_str(r'''
is_irreducible(p) = {
    return polisirreducible(p);
}

get_genus(f) = {
    local(C);
    C = hyperellipticurve(f);
    return hyperellipticgenus(C);
}

get_conductor(f) = {
    local(C);
    C = hyperellipticurve(f);
    return hyperellconductor(C);
}
''')

# Define the polynomial from the curve equation y^2 = f(x)
x = pari.var('x')
f = x**6 + 4*x**5 + 6*x**4 + 2*x**3 + x**2 + 2*x + 1

# Check if the polynomial is irreducible (implies no repeated roots over Q)
is_irred = pari.is_irreducible(f)

# Get the genus of the curve
genus = pari.get_genus(f)

print(f"The polynomial f(x) = {str(f)}.")
if is_irred:
    print("The polynomial is irreducible over the rational numbers.")
else:
    print("The polynomial is reducible, which might lead to a simpler curve.")

print(f"The genus of the curve y^2 = f(x) is {genus}.")

if genus == 1:
    print("The curve is an elliptic curve.")
    # The following would be needed to transform to Weierstrass form and compute conductor
    # This part of the code is for illustration, as the genus is not 1.
    # E = pari.ellinit(f)
    # conductor = E.conductor
    # print(f"The conductor of the elliptic curve is {conductor}.")
else:
    print("The curve is not an elliptic curve.")
    print("Assuming the question asks for the conductor of the hyperelliptic curve (of its Jacobian variety).")
    
    # Compute the conductor of the hyperelliptic curve
    conductor = pari.get_conductor(f)
    print(f"\nThe conductor of the hyperelliptic curve is: {conductor}")
