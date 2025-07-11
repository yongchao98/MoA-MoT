import textwrap

def find_triboluminescent_bacterium():
    """
    Explains and identifies the genus of bacterium known for producing a
    compound with triboluminescent properties.
    """
    explanation = """
    The phenomenon in question is triboluminescence, which is the emission of light from a substance when it is mechanically broken or crushed.

    1. A well-known compound that exhibits triboluminescence when crystallized is monosodium glutamate (MSG).

    2. MSG is the sodium salt of glutamic acid.

    3. The bacterium Corynebacterium glutamicum is used extensively in industrial processes to produce glutamic acid through fermentation.

    4. Therefore, the genus Corynebacterium is known for producing a chemical compound (glutamic acid) that, when converted to a salt and crystallized, exhibits triboluminescent properties.
    """
    
    answer_genus = "Corynebacterium"
    
    print(textwrap.dedent(explanation).strip())
    print("\n---")
    print(f"The genus of bacterium is: {answer_genus}")

if __name__ == "__main__":
    find_triboluminescent_bacterium()