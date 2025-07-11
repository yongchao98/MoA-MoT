import sys

def get_product_name():
  """
  This function returns the IUPAC name of the major product.
  The reaction is the thermal syn-elimination of ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene.
  The sulfoxide elimination yields an alkene and a sulfenic acid.
  Reactant: Ph-S(=O)-CH2-CH2-O-C(CH3)2-CH=CH2
  Product: CH2=CH-O-C(CH3)2-CH=CH2
  The IUPAC name for this product is 3-(ethenyloxy)-3-methylbut-1-ene.
  The numbers in the name are 3, 3, and 1.
  """
  product_name = "3-(ethenyloxy)-3-methylbut-1-ene"
  return product_name

# Print the final IUPAC name.
# The numbers in the name are part of the output as requested.
print(get_product_name())
