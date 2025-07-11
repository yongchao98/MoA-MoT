import sys

# Suppress the creation of __pycache__
sys.dont_write_bytecode = True

class KinshipSystem:
    """A class to represent and compare kinship systems based on Lévi-Strauss's model."""
    def __init__(self, name, lineage, fs_son, mb_ss, hw_wife, bro_sis):
        self.name = name
        self.lineage = lineage  # 'patrilineal' or 'matrilineal'
        self.fs_son = fs_son    # Father-Son relationship (+/-)
        self.mb_ss = mb_ss      # Mother's Brother - Sister's Son (+/-)
        self.hw_wife = hw_wife  # Husband-Wife relationship (+/-)
        self.bro_sis = bro_sis  # Brother-Sister relationship (+/-)

    def matches(self, other):
        """Checks if this system's structure matches another."""
        return (self.lineage == other.lineage and
                self.fs_son == other.fs_son and
                self.mb_ss == other.mb_ss and
                self.hw_wife == other.hw_wife and
                self.bro_sis == other.bro_sis)

def solve_kinship_diagram():
    """Interprets the diagram and evaluates the answer choices."""
    # Step 1: Interpret the Lévi-Strauss diagram based on the image
    # The diagram shows:
    # - A man (Δ), his wife (o), and their son (Δ below).
    # - The wife's brother (Δ to the left).
    # - Attitudes are marked with '+' (familiar) and '-' (formal).
    # - Father-Son relationship is marked with '-': formal/respectful.
    # - Mother's Brother-Sister's Son relationship is marked with '+': familiar/indulgent.
    # - Husband-Wife relationship is marked with '+': familiar.
    # - Brother-Sister relationship is marked with '-': formal/avoidance.
    
    # A system where authority (formal relationship) lies with the father is patrilineal.
    diagram_system = KinshipSystem(
        name="Diagram",
        lineage="patrilineal",
        fs_son="-",
        mb_ss="+",
        hw_wife="+",
        bro_sis="-"
    )

    print("Step 1: Analyzing the Kinship Diagram")
    print("The diagram represents a system with the following characteristics:")
    print(f"- Lineage: {diagram_system.lineage} (Authority lies with the Father, not the Maternal Uncle).")
    print(f"- Father-Son Attitude: Formal/Antagonistic ('{diagram_system.fs_son}')")
    print(f"- Maternal Uncle-Nephew Attitude: Familiar/Tender ('{diagram_system.mb_ss}')")
    print(f"- Husband-Wife Attitude: Familiar/Tender ('{diagram_system.hw_wife}')")
    print(f"- Brother-Sister Attitude: Formal/Avoidance ('{diagram_system.bro_sis}')")
    print("\nWe must find the pair of societies that fits this specific patrilineal structure.")
    print("-" * 20)

    # Step 2: Define the kinship systems from the answer choices based on anthropological data.
    # The 'matrilineal' type is the inverse of the 'patrilineal' on the key axes.
    # The 'Cherkess' system represents a different type of patrilineal structure.
    systems = {
        "Trobriand": KinshipSystem("Trobriand", "matrilineal", fs_son="+", mb_ss="-", hw_wife="-", bro_sis="+"),
        "Siuoi": KinshipSystem("Siuoi", "matrilineal", fs_son="+", mb_ss="-", hw_wife="-", bro_sis="+"),
        "Lake Kubutu": KinshipSystem("Lake Kubutu", "patrilineal", fs_son="-", mb_ss="+", hw_wife="+", bro_sis="-"),
        "Tonga": KinshipSystem("Tonga", "patrilineal", fs_son="-", mb_ss="+", hw_wife="+", bro_sis="-"),
        "Cherkess": KinshipSystem("Cherkess", "patrilineal", fs_son="-", mb_ss="+", hw_wife="-", bro_sis="+")
    }

    choices = {
        "A": ["Trobriand", "Siuoi"],
        "B": ["Siuoi", "Lake Kubutu"],
        "C": ["Lake Kubutu", "Tonga"],
        "D": ["Tonga", "Cherkess"],
        "E": ["Cherkess", "Trobriand"]
    }
    
    print("Step 2: Evaluating Each Answer Choice\n")
    correct_choice_letter = ""

    for choice, societies in choices.items():
        s1_name, s2_name = societies[0], societies[1]
        system1 = systems[s1_name]
        system2 = systems[s2_name]

        match1 = system1.matches(diagram_system)
        match2 = system2.matches(diagram_system)

        print(f"Checking Choice {choice}: {s1_name} and {s2_name}")
        
        # Analyze first society in the pair
        if match1:
            print(f"- {s1_name}: MATCHES the diagram's structure.")
        else:
            print(f"- {s1_name}: DOES NOT MATCH.")
            if system1.lineage != diagram_system.lineage:
                print(f"  Reason: Incorrect lineage type ({system1.lineage}).")
            else:
                print(f"  Reason: Incorrect pattern of intimate relations (H-W is '{system1.hw_wife}', B-S is '{system1.bro_sis}').")
        
        # Analyze second society in the pair
        if match2:
            print(f"- {s2_name}: MATCHES the diagram's structure.")
        else:
            print(f"- {s2_name}: DOES NOT MATCH.")
            if system2.lineage != diagram_system.lineage:
                 print(f"  Reason: Incorrect lineage type ({system2.lineage}).")
            else:
                print(f"  Reason: Incorrect pattern of intimate relations (H-W is '{system2.hw_wife}', B-S is '{system2.bro_sis}').")

        if match1 and match2:
            correct_choice_letter = choice
            print(f"Result: Choice {choice} is the correct answer.\n")
        else:
            print(f"Result: Choice {choice} is incorrect.\n")
        print("-" * 20)

solve_kinship_diagram()