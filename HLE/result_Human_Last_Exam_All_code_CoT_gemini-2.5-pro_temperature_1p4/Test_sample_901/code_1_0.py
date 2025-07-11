# The reaction is an E2 elimination of (1S,2R)-1-bromo-2-methylcyclohexane with potassium tert-butoxide.
# The E2 mechanism requires an anti-periplanar (trans-diaxial) arrangement of the leaving group (Br) and a beta-proton.
# The reaction must proceed through the conformer where Br is axial.
# In the (1S,2R) isomer, when Br is axial, the 2-methyl group is also axial.
# The beta-protons are on C2 and C6.
# - The C2 proton is equatorial (since the methyl group is axial). It cannot be eliminated. This would have formed the Zaitsev product (1-methylcyclohexene).
# - The C6 proton is axial and is anti-periplanar to the axial Br.
# Therefore, elimination occurs by removing the C6 proton to form a double bond between C1 and C6.
# The resulting product is 3-methylcyclohexene.

# The final equation is the name of the product.
# We will construct the name and output it. The number in the name is 3.
product_locant = 3
product_name_base = "methylcyclohexene"
product_name = f"{product_locant}-{product_name_base}"

print("The name of the major product is:")
print(product_name)