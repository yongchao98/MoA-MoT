import sys
import io

# Set stdout to utf-8 to ensure pi symbol prints correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def print_cycloaddition_possibilities():
    """
    Prints four possible [mπ+nπ] descriptions for the dimerization of
    3-oxidopyrylium.
    """
    pi_symbol = "\u03c0"

    # The 3-oxidopyrylium monomer can be formally treated as a 2π, 4π, or 6π
    # component in a cycloaddition. In a dimerization, two such monomers
    # react, leading to several possible formal descriptions.

    # Possibility 1: [4π+2π]
    # One molecule acts as a 4π diene, the other as a 2π dienophile.
    m1, n1 = 4, 2
    print(f"Possibility 1: [{m1}{pi_symbol} + {n1}{pi_symbol}] cycloaddition")

    # Possibility 2: [6π+2π]
    # One molecule acts as a 6π component, the other as a 2π component.
    m2, n2 = 6, 2
    print(f"Possibility 2: [{m2}{pi_symbol} + {n2}{pi_symbol}] cycloaddition")

    # Possibility 3: [4π+4π]
    # Both molecules act as 4π components (e.g., diene + diene).
    m3, n3 = 4, 4
    print(f"Possibility 3: [{m3}{pi_symbol} + {n3}{pi_symbol}] cycloaddition")

    # Possibility 4: [6π+4π]
    # One molecule acts as a 6π component, the other as a 4π component.
    m4, n4 = 6, 4
    print(f"Possibility 4: [{m4}{pi_symbol} + {n4}{pi_symbol}] cycloaddition")

if __name__ == '__main__':
    print_cycloaddition_possibilities()
