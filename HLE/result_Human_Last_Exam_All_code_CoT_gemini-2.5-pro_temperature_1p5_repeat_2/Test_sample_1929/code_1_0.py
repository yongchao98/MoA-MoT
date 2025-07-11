import re

def main():
    """
    This script describes the chemical reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and identifies the product.
    """

    # --- Reactants ---
    diene = {"name": "Butadiene", "formula": "C4H6"}
    dienophile = {"name": "1,1-dichloro-2,2-difluoroethene", "formula": "C2Cl2F2"}

    # --- Product ---
    product = {"name": "4,4-dichloro-5,5-difluorocyclohex-1-ene", "formula": "C6H6Cl2F2"}

    # --- Print Reaction Summary ---
    print("This program identifies the product of a Diels-Alder reaction.\n")
    print("--- Reaction Equation ---")
    print(f"{diene['name']} ({diene['formula']}) + {dienophile['name']} ({dienophile['formula']}) -> {product['name']} ({product['formula']})")
    print("-" * 25)

    print("\n--- Explanation ---")
    print("The reaction is a [4+2] cycloaddition (Diels-Alder reaction) that forms a six-membered ring.")
    print(f"The final product is named: {product['name']}")
    print("-" * 25)

    # --- Output numbers from the final IUPAC name ---
    print("\nThe numbers in the final product's IUPAC name are:")
    # Use regex to find all numbers in the product name
    numbers = re.findall(r'\d+', product['name'])
    for num in numbers:
      print(num)

if __name__ == "__main__":
    main()