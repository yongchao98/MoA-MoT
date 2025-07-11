def solve():
    """
    Provides the scansion for the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    """
    # This line from Virgil's Aeneid (6.881) is metrically complex.
    # The accepted scansion requires applying several special rules like elision,
    # synizesis, and hiatus.
    # A dactyl is represented by "D" (long-short-short).
    # A spondee is represented by "S" (long-long).

    # The consensus scansion for the six feet is:
    # Dactyl, Spondee, Dactyl, Dactyl, Dactyl, Spondee
    scansion_feet = ["D", "S", "D", "D", "D", "S"]

    print("Scansion of 'verum ubi equi atque hominis casu convenit imago.'")
    # Print each foot's representation
    for foot in scansion_feet:
        print(foot, end=" ")
    print() # for a final newline

solve()
<<<D S D D D S>>>